#include <fstream>
#include <vector>
#include <fem/DofSpace.h>

#include <iostream>

DofSpace::DofSpace(std::shared_ptr<ElementSet> elements, std::shared_ptr<NodeSet> nodes)
{
  m_dofTypes = elements->GetDofType();
  m_nodeCoords = &nodes->m_nodeCoords;

  int count = 0, iNode = 0;
  std::vector<int> dofIndex;
  for(auto node : (*m_nodeCoords))
  {
    m_IDmap.emplace_back(node.first);
    dofIndex.clear();
    for(int j = 0; j < m_dofTypes.size(); j++){
      dofIndex.emplace_back(count);
      count += 1;
    }
    m_dofs.emplace_back(dofIndex);
  }
}

DofSpace::~DofSpace()
{}

void DofSpace::ReadFromFile(const std::string&fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  std::string line = "";
  while(true)
  {
    getline(fin, line); line.erase(line.find("\r"));

    if(line.npos != line.find("<NodeConstraints>"))
    {
      while(true)
      {
        getline(fin, line); line.erase(line.find("\r"));
        line.erase(0, line.find_first_not_of(" "));

        if(line.npos != line.find("</NodeConstraints>")) return;

        std::string tempStr = Tools::StringStrip(line);
        std::vector<std::string> a = Tools::StringSplit(tempStr, ";");

        std::vector<std::string> b = Tools::StringSplit(a[0], "=");
        if(2 == b.size()){
          std::vector<std::string> c = Tools::StringSplit(b[0], "[");
          std::string dofType = c[0];
          int nodeID = std::stoi(Tools::StringSplit(c[1],"]")[0]);
          
          Constrain(nodeID, dofType, std::stod(b[1]));
        }
      }
    }
  }
}

void DofSpace::Constrain(const int&nodeId, const std::string&dofType, const double&value)
{
  if(0 == m_nodeCoords->count(nodeId)) {
    throw "Node ID " + std::to_string(nodeId) + " does not exist";
  }

  if(m_dofTypes.end() == std::find(m_dofTypes.begin(), m_dofTypes.end(), dofType)) {
    throw "DOF type " + dofType + " does not exit";
  }
  int indexRow = std::distance(m_IDmap.begin(), std::find(m_IDmap.begin(), m_IDmap.end(), nodeId));
  int indexLine = std::distance(m_dofTypes.begin(), std::find(m_dofTypes.begin(), m_dofTypes.end(), dofType));
  m_constrained.insert(std::pair<int, double>(m_dofs[indexRow][indexLine], value));
}


PetscErrorCode DofSpace::Solve(Mat&K, Vec&df, Vec&da, KSP&ksp)
{
  PetscErrorCode ierr;
  Mat C;
  GetConstraintsMatrix(C);
  // ierr = MatView(C, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  Vec a;
  ierr = VecCreate(PETSC_COMM_WORLD, &a); CHKERRQ(ierr);
  ierr = VecSetSizes(a, PETSC_DECIDE, m_dofs.size()*m_dofs[0].size()); CHKERRQ(ierr);
  ierr = VecSetFromOptions(a); CHKERRQ(ierr);

  for(auto iConstrained : m_constrained)
    ierr = VecSetValue(a, iConstrained.first, m_constrainedFac*iConstrained.second, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(a); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(a); CHKERRQ(ierr);

  Mat constrainedK;
  int numOfRow = 0, numOfLine = 0;
  ierr = MatGetSize(C, &numOfRow, &numOfLine); CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD, &constrainedK); CHKERRQ(ierr);
  ierr = MatSetSizes(constrainedK, PETSC_DECIDE, PETSC_DECIDE, numOfLine, numOfLine); CHKERRQ(ierr);

  Mat tempMat;
  ierr = MatTransposeMatMult(C, K, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tempMat); CHKERRQ(ierr);
  ierr = MatMatMult(tempMat, C, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &constrainedK); CHKERRQ(ierr);
  // std::cout << "constrainedK = " << std::endl;
  // ierr = MatView(constrainedK, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  Vec constrainedB;
  Vec tempVec1, tempVec2;
  ierr = VecCreate(PETSC_COMM_WORLD, &constrainedB); CHKERRQ(ierr);
  ierr = VecSetSizes(constrainedB, PETSC_DECIDE, numOfLine); CHKERRQ(ierr);
  ierr = VecSetFromOptions(constrainedB); CHKERRQ(ierr);
  
  ierr = VecDuplicate(a, &tempVec1); CHKERRQ(ierr);
  ierr = VecAXPY(tempVec1, -1., a); CHKERRQ(ierr);
  ierr = VecDuplicate(a, &tempVec2); CHKERRQ(ierr);
  // ierr = VecScale(tempVec1, -1.); CHKERRQ(ierr);
  // std::cout << "tempVec1 = " << std::endl;
  // ierr = VecView(tempVec1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = MatMult(K, tempVec1, tempVec2); CHKERRQ(ierr);
  // std::cout << "tempVec2 = " << std::endl;
  // ierr = VecView(tempVec2, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // std::cout << "df = " << std::endl;
  // ierr = VecView(df, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecWAXPY(tempVec1, 1., df, tempVec2); CHKERRQ(ierr);
  // ierr = VecView(tempVec1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = MatMultTranspose(C, tempVec1, constrainedB); CHKERRQ(ierr);
  // std::cout << "constrainedB = " << std::endl;
  // ierr = VecView(constrainedB, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  Vec constrainedDa;
  ierr = VecDuplicate(constrainedB, &constrainedDa); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, constrainedK, constrainedK); CHKERRQ(ierr);
  ierr = KSPSetUp(ksp); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, constrainedB, constrainedDa); CHKERRQ(ierr);
  // std::cout << "constrainedDa = " << std::endl;
  // ierr = VecView(constrainedDa, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = MatMult(C, constrainedDa, da); CHKERRQ(ierr);
  ierr = VecAXPY(da, 1.0, a); CHKERRQ(ierr);

  ierr = MatDestroy(&C); CHKERRQ(ierr); 
  ierr = MatDestroy(&constrainedK); CHKERRQ(ierr);
  ierr = MatDestroy(&tempMat); CHKERRQ(ierr);

  ierr = VecDestroy(&tempVec1); CHKERRQ(ierr);
  ierr = VecDestroy(&tempVec2); CHKERRQ(ierr);
  ierr = VecDestroy(&constrainedDa); CHKERRQ(ierr);
  ierr = VecDestroy(&constrainedB); CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode DofSpace::Solve(double K, Vec&df, Vec&da)
{
  PetscErrorCode ierr;
  if (abs(K) > 1e-10)
  {
    Vec temp;
    ierr = VecDuplicate(da, &temp); CHKERRQ(ierr);
    ierr = VecSet(temp, 0.); CHKERRQ(ierr);
    ierr = VecWAXPY(da, 1. / K, df, temp); CHKERRQ(ierr);
    ierr = VecDestroy(&temp); CHKERRQ(ierr);
  }
  else
  {
    throw "denominator is too small";
  }
  return ierr;
}

PetscErrorCode DofSpace::GetConstraintsMatrix(Mat &C)
{
  int n_constrianed = m_constrained.size();
  int n = m_dofs.size() * m_dofs[0].size();

  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, &C); CHKERRQ(ierr);
  ierr = MatSetSizes(C, PETSC_DECIDE, PETSC_DECIDE, n, n-n_constrianed); CHKERRQ(ierr);
  ierr = MatSetUp(C); CHKERRQ(ierr);
  
  int j = 0;
  for(int i = 0; i < n ; i++)
  {
    if(0 != m_constrained.count(i)) continue;
    ierr = MatSetValue(C, i, j, 1., INSERT_VALUES); CHKERRQ(ierr);
    j += 1;
  }
  ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode DofSpace::Norm(Vec &r, double &error)
{
  PetscErrorCode ierr;

  Mat C;
  GetConstraintsMatrix(C);

  int numOfLine;
  ierr = MatGetSize(C, nullptr, &numOfLine); CHKERRQ(ierr);

  Vec tempVec;
  ierr = VecCreate(PETSC_COMM_WORLD, &tempVec); CHKERRQ(ierr);
  ierr = VecSetSizes(tempVec, PETSC_DECIDE, numOfLine); CHKERRQ(ierr);
  ierr = VecSetFromOptions(tempVec); CHKERRQ(ierr);

  ierr = MatMultTranspose(C, r, tempVec); CHKERRQ(ierr);
  ierr = VecNorm(tempVec, NORM_2, &error); CHKERRQ(ierr);

  ierr = VecDestroy(&tempVec);

  return ierr;
}