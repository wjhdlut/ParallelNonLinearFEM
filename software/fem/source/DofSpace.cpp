/**
 * @File Name:     DofSpace.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <regex>

#include "../include/DofSpace.h"
#include "../../util/include/DataStructure.h"


DofSpace::DofSpace(std::shared_ptr<ElementSet> elements, std::shared_ptr<NodeSet> nodes)
{
  m_dofTypes = elements->GetDofType();
  m_nodeCoords = &nodes->m_nodeCoords;

  Initialize();
}

DofSpace::~DofSpace()
{
  if(NULL != m_transMatrix) MatDestroy(&m_transMatrix);
  if(NULL != m_C) MatDestroy(&m_C);
}

void DofSpace::Initialize()
{
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

void DofSpace::ReadFromFile(const std::string&fileName)
{
  ReadNodeConstraint(fileName);

  ReadRigidWall(fileName);

  ReadMasterSlaveConstraint(fileName);
}

void DofSpace::Constrain(const int&nodeId, const std::string&dofType, const double&value)
{
  if(0 == m_nodeCoords->count(nodeId)) {
    std::cout << "Catch Exception: "
              << "Node ID " + std::to_string(nodeId) + " does not exist"
              << std::endl;
    exit(-1);
  }

  if(m_dofTypes.end() == std::find(m_dofTypes.begin(), m_dofTypes.end(), dofType)) {
    std::cout << "Catch Exception: "
              << "DOF type " + dofType + " does not exit"
              << std::endl;
    exit(-1);
  }
  int indexRow = std::distance(m_IDmap.begin(), std::find(m_IDmap.begin(), m_IDmap.end(), nodeId));
  int indexLine = std::distance(m_dofTypes.begin(), std::find(m_dofTypes.begin(), m_dofTypes.end(), dofType));
  m_constrained.insert(std::pair<int, double>(m_dofs[indexRow][indexLine], value));
}

PetscErrorCode DofSpace::Solve(Mat&K, Vec&df, Vec&da, KSP&ksp)
{
  // Transform the Stiffness Matrix and Force Vector 
  // from Cartesian Coordinate System to Cylindrical Coordinate System 
  if(m_transMatrix != NULL)
  {
    TransMatToCylinder(K);
    TransVecToCylinder(df);
  }
  
  PetscErrorCode ierr;
  // Mat C;
  // GetConstraintsMatrix(C);
  // ierr = MatView(C, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  Vec a;
  ierr = VecCreate(PETSC_COMM_WORLD, &a); CHKERRQ(ierr);
  ierr = VecSetSizes(a, PETSC_DECIDE, m_dofs.size()*m_dofs[0].size()); CHKERRQ(ierr);
  ierr = VecSetFromOptions(a); CHKERRQ(ierr);

  m_constrainedFac = GlobalData::GetInstance()->m_lam;
  for(auto iConstrained : m_constrained)
    ierr = VecSetValue(a, iConstrained.first, m_constrainedFac*iConstrained.second, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(a); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(a); CHKERRQ(ierr);
  std::string fileA="a.txt";
  Tools::PrintVecIntoFile(a, fileA);

  Mat constrainedK;
  int numOfRow = 0, numOfLine = 0;
  ierr = MatGetSize(m_C, &numOfRow, &numOfLine); CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD, &constrainedK); CHKERRQ(ierr);
  ierr = MatSetSizes(constrainedK, PETSC_DECIDE, PETSC_DECIDE, numOfLine, numOfLine); CHKERRQ(ierr);

  Mat tempMat;
  ierr = MatTransposeMatMult(m_C, K, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tempMat); CHKERRQ(ierr);
  ierr = MatMatMult(tempMat, m_C, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &constrainedK); CHKERRQ(ierr);
  // std::cout << "constrainedK = " << std::endl;
  // ierr = MatView(constrainedK, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // std::string fileK = "constrainedK.txt";
  // ierr = Tools::PrintMatIntoFile(constrainedK, fileK);

  Vec constrainedB;
  Vec tempVec1, tempVec2;
  ierr = VecCreate(PETSC_COMM_WORLD, &constrainedB); CHKERRQ(ierr);
  ierr = VecSetSizes(constrainedB, PETSC_DECIDE, numOfLine); CHKERRQ(ierr);
  ierr = VecSetFromOptions(constrainedB); CHKERRQ(ierr);
  
  ierr = VecDuplicate(a, &tempVec1); CHKERRQ(ierr);
  ierr = VecZeroEntries(tempVec1); CHKERRQ(ierr);
  ierr = VecAXPY(tempVec1, -1., a); CHKERRQ(ierr);
  ierr = VecDuplicate(a, &tempVec2); CHKERRQ(ierr);
  ierr = VecZeroEntries(tempVec2); CHKERRQ(ierr);
  // ierr = VecScale(tempVec1, -1.); CHKERRQ(ierr);
  // std::cout << "tempVec1 = " << std::endl;
  // ierr = VecView(tempVec1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // f = K*u;
  ierr = MatMult(K, tempVec1, tempVec2); CHKERRQ(ierr);
  // std::cout << "tempVec2 = " << std::endl;
  // ierr = VecView(tempVec2, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  // std::cout << "df = " << std::endl;
  // ierr = VecView(df, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecWAXPY(tempVec1, 1., df, tempVec2); CHKERRQ(ierr);
  // ierr = VecView(tempVec1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  
  ierr = MatMultTranspose(m_C, tempVec1, constrainedB); CHKERRQ(ierr);
  // std::cout << "constrainedB = " << std::endl;
  // ierr = VecView(constrainedB, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  std::string fileB = "constrainedB.txt";
  ierr = Tools::PrintVecIntoFile(constrainedB, fileB);

  Vec constrainedDa;
  ierr = VecDuplicate(constrainedB, &constrainedDa); CHKERRQ(ierr);

  PC kspPC;
  ierr = KSPSetType(ksp, KSPFCG); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &kspPC); CHKERRQ(ierr);
  ierr = PCSetType(kspPC, PCSOR); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, constrainedK, constrainedK); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSetUp(ksp); CHKERRQ(ierr);
  // ierr = KSPSetType(ksp, KSPFCG); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, constrainedB, constrainedDa); CHKERRQ(ierr);
  // std::cout << "constrainedDa = " << std::endl;
  // ierr = VecView(constrainedDa, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // Tools::PrintVecIntoFile(constrainedDa, "constrainedDa.txt");
  ierr = MatMult(m_C, constrainedDa, da); CHKERRQ(ierr);
  ierr = VecAXPY(da, 1.0, a); CHKERRQ(ierr);
  // std::cout << "da = " << std::endl;
  // ierr = VecView(da, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // Tools::PrintVecIntoFile(da, "da.txt");
  
  if(m_transMatrix != NULL){
    ierr = TransVecToGlobalCSY(da); CHKERRQ(ierr);
  }

  ierr = MatDestroy(&constrainedK); CHKERRQ(ierr);
  ierr = MatDestroy(&tempMat); CHKERRQ(ierr);

  ierr = VecDestroy(&tempVec1); CHKERRQ(ierr);
  ierr = VecDestroy(&tempVec2); CHKERRQ(ierr);
  ierr = VecDestroy(&constrainedDa); CHKERRQ(ierr);
  ierr = VecDestroy(&constrainedB); CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode DofSpace::Solve(const Vec &K, const Vec &df, Vec&da)
{
  PetscErrorCode ierr;
  
  // computed acce
  ierr = VecPointwiseDivide(da, df, K); CHKERRQ(ierr);
  for(auto iConstrained : m_constrained)
    ierr = VecSetValue(da, iConstrained.first, m_constrainedFac*iConstrained.second, INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(da); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(da); CHKERRQ(ierr);
  // Tools::PrintVecIntoFile(da, "da.txt");

  if(nullptr != m_rigidWall) RigidWallConstraint(da);

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

PetscErrorCode DofSpace::GetConstraintsMatrix()
{
  int n_constrianed = m_constrained.size();
  int n = m_dofs.size() * m_dofs[0].size();
  int numOfSlave = 0;
  for(auto iMasterSlave : m_masterSlave){
    for(auto iSlave : iMasterSlave.second){
      if(0 != m_constrained.count(iSlave)) continue;
      numOfSlave += 1;
    }
  }

  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, &m_C); CHKERRQ(ierr);
  ierr = MatSetSizes(m_C, PETSC_DECIDE, PETSC_DECIDE, n, n-n_constrianed-numOfSlave); CHKERRQ(ierr);
  ierr = MatSetUp(m_C); CHKERRQ(ierr);
  
  int j = 0;
  for(int i = 0; i < n ; i++)
  {
    bool pass = false;
    for(auto iMasterSlave : m_masterSlave){
      if(0 != iMasterSlave.second.count(i)){
        pass = true;
        ierr = MatSetValue(m_C, i, iMasterSlave.first, 1., INSERT_VALUES); CHKERRQ(ierr);
        break;
      }
    }
    if(0 != m_constrained.count(i) || pass) continue;
    ierr = MatSetValue(m_C, i, j, 1., INSERT_VALUES); CHKERRQ(ierr);
    j += 1;
  }
  ierr = MatAssemblyBegin(m_C, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_C, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  Tools::PrintMatIntoFile(m_C, "transMatrix.txt");
  
  return ierr;
}

PetscErrorCode DofSpace::Norm(const Vec &r, double &error)
{
  PetscErrorCode ierr;

  // Mat C;
  // GetConstraintsMatrix(C);
  Vec rCopy;
  ierr = VecDuplicate(r, &rCopy); CHKERRQ(ierr);
  ierr = VecCopy(r, rCopy); CHKERRQ(ierr);
  if(m_transMatrix != NULL)
    TransVecToCylinder(rCopy);

  int numOfLine;
  ierr = MatGetSize(m_C, nullptr, &numOfLine); CHKERRQ(ierr);

  Vec tempVec;
  ierr = VecCreate(PETSC_COMM_WORLD, &tempVec); CHKERRQ(ierr);
  ierr = VecSetSizes(tempVec, PETSC_DECIDE, numOfLine); CHKERRQ(ierr);
  ierr = VecSetFromOptions(tempVec); CHKERRQ(ierr);

  ierr = MatMultTranspose(m_C, rCopy, tempVec); CHKERRQ(ierr);
  ierr = VecNorm(tempVec, NORM_2, &error); CHKERRQ(ierr);

  // std::cout << "tempVec = \n" << std::endl;
  // ierr = VecView(tempVec, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = VecDestroy(&tempVec); CHKERRQ(ierr);
  ierr = VecDestroy(&rCopy); CHKERRQ(ierr);

  return ierr;
}

void DofSpace::ReadNodeConstraint(const std::string &fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  std::string line = "";
  while(true)
  {
    getline(fin, line);
    if(line.npos != line.find("\r"))
      line.erase(line.find("\r"));

    if(line.npos != line.find("<NodeConstraints>"))
    {
      while(true)
      {
        getline(fin, line);
        if(line.npos !=line.find("\r"))
          line.erase(line.find("\r"));
        line.erase(0, line.find_first_not_of(" "));
        
        if(line.size() == 0) continue;

        if(line.npos != line.find("</NodeConstraints>")){
          fin.close();
          return;
        }

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
    if(fin.eof()) break;
  }
  fin.close();
}

void DofSpace::ReadRigidWall(const std::string &fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  std::string line = "";
  while(true)
  {
    getline(fin, line);
    if(line.npos != line.find("\r"))
      line.erase(line.find("\r"));

    if(line.npos != line.find("<RigidWall>"))
    {
      while(true)
      {
        getline(fin, line);
        if(line.npos !=line.find("\r"))
          line.erase(line.find("\r"));
        line.erase(0, line.find_first_not_of(" "));
        
        if(line.size() == 0) continue;

        if(line.npos != line.find("</RigidWall>")) {
          fin.close();
          return;
        }

        std::string tempStr = Tools::StringStrip(line);
        std::vector<std::string> a = Tools::StringSplit(tempStr, ";");

        std::vector<std::string> b = Tools::StringSplit(a[0], "=");
        
        if(2 > b.size()){
          std::cout << "Catch Exception: "
                    << "The input of RigidWall is wrong!"
                    << std::endl;
          exit(-1);
        }
        m_rigidWall = std::make_shared<RigidWall>();
        if(std::string::npos != b[0].find('x')) m_rigidWall->direction = 1;
        if(std::string::npos != b[0].find('y')) m_rigidWall->direction = 2;
        if(std::string::npos != b[0].find('z')) m_rigidWall->direction = 3;
        m_rigidWall->coord = std::stod(b[1]);
      }
    }
    if(fin.eof()) break;
  }
  fin.close();
}

void DofSpace::ReadMasterSlaveConstraint(const std::string &fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  std::string line = "";
  std::regex pattern("[\\s]{2,}");
  while(true)
  {
    getline(fin, line);
    if(line.npos != line.find("\r"))
      line.erase(line.find("\r"));

    if(fin.eof()) break;

    if(line.npos != line.find("<MasterSlaveNodalConstraints>"))
    {
      while(true)
      {
        getline(fin, line);
        if(line.npos != line.find("\r"))
          line.erase(line.find("\r"));
        line.erase(0, line.find_first_not_of(" "));

        if(line.npos != line.find("</MasterSlaveNodalConstraints>"))
        {
          fin.close();
          return;
        }
        
        line =  std::regex_replace(line, pattern, " ");
        std::vector<std::string> a = Tools::StringSplit(line, ";");
        std::vector<std::string> b = Tools::StringSplit(a[0], " ");

        if(3 > b.size()){
          std::cout << "Catch Exception: "
                    << "The input of MasterSlaveConstraint is wrong!"
                    << std::endl;
          exit(-1);
        }

        int masterNodeId = std::stoi(b[0]);
        std::set<int> slaveNodesId;
        for(int i = 2; i < b.size(); i++)
          slaveNodesId.insert(std::stoi(b[i]));
        MasterSlaveConstraint(masterNodeId, b[1], slaveNodesId);
      }
    }
  }
  fin.close();
}

void DofSpace::MasterSlaveConstraint(const int &masterNodeId,
                                     const std::string &dofType,
                                     const std::set<int> &slaveNodesId)
{
  if(0 == m_nodeCoords->count(masterNodeId)){
    std::cout << "Catch Exception: "
              << "Node ID " + std::to_string(masterNodeId) + " does not exist"
              << std::endl;
    exit(-1);
  }
  
  for(int i = 1; i < dofType.size()-1; i++)
  {
    if(m_dofTypes.end() == std::find(m_dofTypes.begin(), m_dofTypes.end(), std::string(1, dofType[i]))){
      std::cout << "Catch Exception: "
                << "DOF type " + dofType + " does not exit"
                << std::endl;
      exit(-1);
    }
    int masterNodeRow  = std::distance(m_IDmap.begin(),
                         std::find(m_IDmap.begin(), m_IDmap.end(), masterNodeId));
    int masterNodeLine = std::distance(m_dofTypes.begin(),
                         std::find(m_dofTypes.begin(), m_dofTypes.end(), std::string(1, dofType[i])));
    int masterDof = m_dofs[masterNodeRow][masterNodeLine];
    
    if(0 != m_masterSlave.count(masterDof)){
      std::cout << " PLEASE CHECK THE MASTER/SLAVE CONSTRAIN FOR MATER NODE ID IS " 
                << masterNodeId
                << std::endl;
    }
    m_masterSlave.insert(std::pair<int, std::set<int>>(masterDof, {}));
    for(auto iter : slaveNodesId)
    {
      int slaveNodeRow = std::distance(m_IDmap.begin(), std::find(m_IDmap.begin(), m_IDmap.end(), iter));
      int slaveNodeDof = m_dofs[slaveNodeRow][masterNodeLine];
      m_masterSlave[masterDof].insert(slaveNodeDof);
    }
  }
}

void DofSpace::RigidWallConstraint(Vec&da)
{
  if(0 != m_rigidWall->direction)
  {
    int index = std::abs(m_rigidWall->direction) - 1;
    double temp = (m_rigidWall->direction < 0) ? -1. : 1.;
    Vec &velo = GlobalData::GetInstance()->m_velo;
    std::map<int, VectorXd> &nodeCoords = GlobalData::GetInstance()->m_nodes->m_nodeCoords;
    const Vec &disp = GlobalData::GetInstance()->m_state;
    
    std::vector<int> dofIndex;
    VectorXd iNodeDisp = VectorXd::Zero(m_dofTypes.size());
    std::vector<double> iNodeVelo(m_dofTypes.size(), 0.);
    std::vector<double> iNodeAcce(m_dofTypes.size(), 0.);
    for(auto iNode : nodeCoords){
      dofIndex = m_dofs[GetIndex(iNode.first)];
      VecGetValues(velo, dofIndex.size(), &dofIndex[0], &iNodeVelo[0]);
      VecGetValues(da, dofIndex.size(), &dofIndex[0], &iNodeAcce[0]);
      VecGetValues(disp, dofIndex.size(), &dofIndex[0], &iNodeDisp[0]);
      VectorXd nodeCoordCurrent = iNode.second + iNodeDisp;
      // std::cout << "nodeCoordCurrent = " << nodeCoordCurrent << std::endl;
      // std::cout << temp*(m_rigidWall->coord-nodeCoordCurrent(index)) << std::endl;
      // std::cout << (temp*(m_rigidWall->coord-nodeCoordCurrent(index)) >= 1.e-10) << std::endl;

      if(temp*(m_rigidWall->coord-nodeCoordCurrent(index)) >= 0.){
        if(temp*iNodeVelo[index] <= 0.) iNodeVelo[index] = 0.;
        if(temp*iNodeAcce[index] <= 0.) iNodeAcce[index] = 0.;

        VecSetValues(velo, dofIndex.size(), &dofIndex[0], &iNodeVelo[0], INSERT_VALUES);
        VecSetValues(da, dofIndex.size(), &dofIndex[0], &iNodeAcce[0], INSERT_VALUES);
      }
      VecAssemblyBegin(velo); VecAssemblyEnd(velo);
      VecAssemblyBegin(da); VecAssemblyEnd(da);
    }
  }
}

PetscErrorCode DofSpace::TransMatToCylinder(Mat&A)
{
  PetscErrorCode ierr;

  ierr = MatTransposeMatMult(m_transMatrix, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A); CHKERRQ(ierr);
  ierr = MatMatMult(A, m_transMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A); CHKERRQ(ierr);
  
  return ierr;
}

PetscErrorCode DofSpace::TransVecToCylinder(Vec&B)
{
  PetscErrorCode ierr;
  
  Vec tempB;
  ierr = VecDuplicate(B, &tempB); CHKERRQ(ierr);
  ierr = MatMultTranspose(m_transMatrix, B, tempB); CHKERRQ(ierr);
  ierr = VecCopy(tempB, B); CHKERRQ(ierr);
  ierr = VecDestroy(&tempB); CHKERRQ(ierr);
  
  return ierr;
}

PetscErrorCode DofSpace::TransVecToGlobalCSY(Vec&B)
{
  PetscErrorCode ierr;
  
  Vec tempB;
  ierr = VecDuplicate(B, &tempB); CHKERRQ(ierr);
  ierr = MatMult(m_transMatrix, B, tempB); CHKERRQ(ierr);
  ierr = VecCopy(tempB, B); CHKERRQ(ierr);
  ierr = VecDestroy(&tempB); CHKERRQ(ierr);
  
  return ierr;
}

PetscErrorCode DofSpace::ComputeTransMatrix()
{
  int numOfNode   = m_dofs.size();
  int numOfDof    = m_dofs.at(0).size();
  int numOfTolDof = numOfNode * numOfDof;

  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, &m_transMatrix); CHKERRQ(ierr);
  ierr = MatSetSizes(m_transMatrix, PETSC_DECIDE, PETSC_DECIDE, numOfTolDof, numOfTolDof); CHKERRQ(ierr);
  ierr = MatSetFromOptions(m_transMatrix); CHKERRQ(ierr);
  ierr = MatSetUp(m_transMatrix); CHKERRQ(ierr);
  
  for(int i = 0; i < m_dofs.size(); i++)
  {
    for(int j = 0; j < m_dofs.at(i).size(); j++){
      int index = m_dofs.at(i).at(j);
      
      if(0 != m_constrained.count(index)){
        double theta = GlobalData::GetInstance()->m_nodes->GetNodeAngle(m_IDmap[i]);
        Vector2i rowIndex = {m_dofs.at(i).at(0), m_dofs.at(i).at(1)};
        Vector2i colIndex = {m_dofs.at(i).at(0), m_dofs.at(i).at(1)};
        Vector4d value = {cos(theta), -sin(theta), sin(theta), cos(theta)};
        MatSetValues(m_transMatrix, 2, &rowIndex[0], 2, &colIndex[0], &value[0], INSERT_VALUES);
        j = 2;
        continue;
      }
      MatSetValue(m_transMatrix, index, index, 1.0, INSERT_VALUES);
      }
  }

  ierr = MatAssemblyBegin(m_transMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_transMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  // std::cout << "transform matrix = \n" << std::endl;
  // ierr = MatView(m_transMatrix, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  return ierr;
}

void DofSpace::CompBasicVariable()
{
  if("Cylindrical" == GlobalData::GetInstance()->m_nodes->GetCoorSysType())
    ComputeTransMatrix();
  
  GetConstraintsMatrix();
}

double DofSpace::Converg(const Vec&fext, const Vec &fint)
{
  Vec fextCopy, fintCopy;
  VecDuplicate(fext, &fextCopy); VecCopy(fext, fextCopy);
  VecDuplicate(fint, &fintCopy); VecCopy(fint, fintCopy);

  // std::cout << "fext = " << std::endl;
  // VecView(fextCopy, PETSC_VIEWER_STDOUT_WORLD);
  // std::cout << "fint = " << std::endl;
  // VecView(fintCopy, PETSC_VIEWER_STDOUT_WORLD);

  TransVecToCylinder(fextCopy);
  TransVecToCylinder(fintCopy);
  
  double tempValue = 0.;
  for(auto iConstrain : m_constrained)
  {
    VecGetValues(fintCopy, 1, &iConstrain.first, &tempValue);
    VecSetValue(fextCopy, iConstrain.first, tempValue, INSERT_VALUES);
  }
  VecAssemblyBegin(fextCopy);
  VecAssemblyEnd(fextCopy);

  TransVecToGlobalCSY(fextCopy);
  TransVecToGlobalCSY(fintCopy);

  // std::cout << "fext after Transform" << std::endl;
  // VecView(fextCopy, PETSC_VIEWER_STDOUT_WORLD);
  // std::cout << "fint after Transform" << std::endl;
  // VecView(fintCopy, PETSC_VIEWER_STDOUT_WORLD);

  VecAYPX(fintCopy, -1., fextCopy);

  double normFext, normDf;
  VecNorm(fextCopy, NORM_2, &normFext);
  VecNorm(fintCopy, NORM_2, &normDf);


  VecDestroy(&fextCopy);
  VecDestroy(&fintCopy);

  return normDf / normFext;
}

// PetscErrorCode DofSpace::SetMasterSlaveConstraint(Mat &K, Vec &dF)
// {
//   PetscErrorCode ierr;
  
//   Mat tempMat;
//   ierr = MatTransposeMatMult(m_masterSlaveTrans, K, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tempMat); CHKERRQ(ierr);
//   ierr = MatMatMult(tempMat, m_masterSlaveTrans, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &K); CHKERRQ(ierr);

//   ierr = MatMultTranspose(m_masterSlaveTrans, dF, dF); CHKERRQ(ierr);

//   return ierr;
// }

// PetscErrorCode DofSpace::GetMasterSlaveMatrix()
// {
//   int numOfDof = m_dofs.size() * m_dofs[0].size();
//   int numOfSlave = 0;
//   for(auto iMasterSlave : m_masterSlave){
//     numOfSlave += iMasterSlave.second.size();
//   }

//   PetscErrorCode ierr;
//   ierr = MatCreate(PETSC_COMM_WORLD, &m_masterSlaveTrans); CHKERRQ(ierr);
//   ierr = MatSetSizes(m_masterSlaveTrans, PETSC_DECIDE, PETSC_DECIDE, numOfDof, numOfDof-numOfSlave); CHKERRQ(ierr);
//   ierr = MatSetUp(m_masterSlaveTrans); CHKERRQ(ierr);
  
//   int j = 0;
//   for(int i = 0; i < numOfDof; i++){
//     bool pass = false;
//     for(auto iMasterSlave : m_masterSlave){
//       if(iMasterSlave.second.count(i)){
//         pass = true;
//         ierr = MatSetValue(m_masterSlaveTrans, i, iMasterSlave.first, 1., INSERT_VALUES); CHKERRQ(ierr);
//         break;
//       }
//     }

//     if(pass) continue;
//     ierr = MatSetValue(m_masterSlaveTrans, i, j, 1., INSERT_VALUES); CHKERRQ(ierr);
//     j += 1;
// }

//   ierr = MatAssemblyBegin(m_masterSlaveTrans, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//   ierr = MatAssemblyEnd(m_masterSlaveTrans, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

//   return ierr;
// }