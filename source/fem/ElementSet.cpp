#include <fstream>
#include <regex>
#include <algorithm>
#include <fem/ElementSet.h>
#include <elements/Element.h>
#include <util/ObjectFactory.h>
#include <util/DataStructure.h>

#include <iostream>

ElementSet::ElementSet(std::shared_ptr<NodeSet> &nodes, const nlohmann::json &props)
{
  m_nodes = nodes;
  m_props = props;
}

ElementSet::~ElementSet()
{}

void ElementSet::ReadFromFile(const std::string &fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  if(!fin.is_open()){
    std::cout << fileName << " open failed!!" << std::endl;
    exit(-1);
  }

  std::string line = "";
  std::regex pattern("[\\s]{2,}");
  while(true)
  {
    getline(fin, line);
    if(line.npos != line.find("\r"))
      line.erase(line.find("\r"));

    if(line.npos != line.find("<Elements>"))
    {
      while(true)
      {
        getline(fin, line);
        if(line.npos != line.find("\r"))
          line.erase(line.find("\r"));

        if(line.npos != line.find("</Elements>")) return;

        if(line.size() == 0) continue;
        
        line =  std::regex_replace(line, pattern, " ");
        std::vector<std::string> strData = Tools::StringSplit(line, ";");
        for(auto iter = strData.begin(); iter != strData.end() - 1; iter++)
        {
          std::string tempStr = Tools::StringStrip(*iter);
          std::vector<std::string> b = Tools::StringSplit(tempStr, " ");
          if("//" == b[0].substr(0, 2) || "#" == b[0].substr(0, 1)) break;

          int elemID = std::stoi(b[0]);
          std::vector<int> elementNodes;
          for(auto iterb = b.begin() + 2; iterb != b.end(); iterb++)
            elementNodes.emplace_back(std::stoi(*iterb));
          
          std::string&modelName = b[1];
          modelName.erase(0, 1); modelName.erase(modelName.size() - 1, 1);
          Add(elemID, modelName, elementNodes);
        }
      }
    }
  }

  fin.close();
}

void ElementSet::Add(const int elemId, const std::string &modelName, const std::vector<int> &elementNodes)
{
  if(!m_props.contains(modelName)) throw "Missing properties for model " + modelName;

  nlohmann::json&modelProps = m_props.at(modelName);

  if(!modelProps.contains("type")) throw "Missing type for model " + modelName;

  std::string modelType = modelProps.at("type");

  std::shared_ptr<Element> elem = ObjectFactory::CreateObject<Element>(modelType, elementNodes, modelProps);

  m_nodes->GetNodeCoords(elem->GetNodes());

  m_elem[elemId] = elem;
  if (0 == m_groups.count(modelName)) m_groups[modelName] = {};
  m_groups[modelName].emplace_back(elemId);
}

std::vector<std::string> ElementSet::GetDofType()
{
  std::vector<std::string> dofTypes = {};

  for(auto element : m_elem)
  {
    for(auto dofType : element.second->GetDofType())
      if(std::find(dofTypes.begin(), dofTypes.end(), dofType) == dofTypes.end())
        dofTypes.emplace_back(dofType);
  }

  return dofTypes;
}

PetscErrorCode ElementSet::AssembleMatrix(Mat&A, Vec&B, const int rank, const std::string&action)
{
  int numOfNode = GlobalData::GetInstance()->m_dofs->m_dofs.size();
  int numOfDof = GlobalData::GetInstance()->m_dofs->m_dofs.at(0).size();
  int numOfTolDof = numOfNode * numOfDof;

  // Create PETSc object Mat
  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, numOfTolDof, numOfTolDof); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatSetUp(A); CHKERRQ(ierr);

  // ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // Create PETSc object Vec
  ierr = VecSet(B, 0.);
  // ierr = VecCreate(PETSC_COMM_WORLD, &B); CHKERRQ(ierr);
  // ierr = VecSetSizes(B, PETSC_DECIDE, numOfTolDof); CHKERRQ(ierr);
  // ierr = VecSetFromOptions(B); CHKERRQ(ierr);

  GlobalData::GetInstance()->ResetNodalOutput();

  m_props = GlobalData::GetInstance()->m_props;
  std::shared_ptr<DofSpace> dofs = GlobalData::GetInstance()->m_dofs;
  
  Vec &state = GlobalData::GetInstance()->m_state;
  Vec &Dstate = GlobalData::GetInstance()->m_Dstate;
  Vec &velo   = GlobalData::GetInstance()->m_velo;

  // ierr = VecView(state, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  std::vector<double> elemMatrixValues;
  for(auto elementGroup : m_groups)
  {
    nlohmann::json &elemPrpos = m_props.at(elementGroup.first);
    for(auto element : elementGroup.second)
    {
      elemPtr = m_elem[element];

      // Get the element nodes
      elemNodes = elemPtr->GetNodes();

      // Get the element coordinates
      elemCoords = m_nodes->GetNodeCoords(elemNodes);

      // Get the element degrees of freedom
      elemDofs = dofs->Get(elemNodes);

      // Get the element state
      elemState.resize(elemDofs.size(), 0.), elemDstate.resize(elemDofs.size(), 0.);
      ierr = VecGetValues(state, elemDofs.size(), &elemDofs[0], &elemState[0]); CHKERRQ(ierr);
      ierr = VecGetValues(Dstate, elemDofs.size(), &elemDofs[0], &elemDstate[0]); CHKERRQ(ierr);

      elemData = std::make_shared<ElementData>(elemState, elemDstate);
      elemData->m_coords = elemCoords;
      ierr = VecGetValues(velo, elemDofs.size(), &elemDofs[0], &elemData->m_velo[0]);
      
      elemPtr->MatReset();
      elemPtr->SetTimeIncrementPara(m_dtK1, m_elemDistortion);
      elemPtr->GetTangentStiffness(elemData);
      elemPtr->ReturnTimeIncrementPara(m_dtK1, m_elemDistortion);
      // std::cout << "element = " << element << "m_dtK1 = " << m_dtK1 << " m_elemDistortion = " << m_elemDistortion << std::endl;
      
      for(auto label : elemData->m_outLabel)
        elemPtr->AppendNodalOutput(label, elemData->m_outputData);

      if(1 == rank)
      {
        ierr = VecSetValues(B, elemDofs.size(), &elemDofs[0], &(elemData->m_fint[0]), ADD_VALUES); CHKERRQ(ierr);
      }
      else if(2 == rank && "getTangentStiffness" == action)
      {
        // Assemble Global Stiffness Matrix and Internal Force Vector
        elemMatrixValues = Math::ConvertMatrixToVec(elemData->m_stiff);
        ierr = MatSetValues(A, elemDofs.size(), &elemDofs[0], elemDofs.size(), 
                           &elemDofs[0], &elemMatrixValues[0], ADD_VALUES); CHKERRQ(ierr);        
        ierr = VecSetValues(B, elemDofs.size(), &elemDofs[0], &(elemData->m_fint[0]), ADD_VALUES); CHKERRQ(ierr);
      }
      else if(2 == rank && "getMassMatrix" == action)
      {
        elemPtr->GetMassMatrix(elemData);
        elemMatrixValues = Math::ConvertMatrixToVec(elemData->m_mass);
        ierr = MatSetValues(A, elemDofs.size(), &elemDofs[0], elemDofs.size(),
                           &elemDofs[0], &elemMatrixValues[0], ADD_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(B, elemDofs.size(), &elemDofs[0], &(elemData->m_lumped[0]), ADD_VALUES); CHKERRQ(ierr);
      }
      else
      {
        throw "assemleArray is only implemented for vectors and matrices.";
      }
    }
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  // ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(B); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(B); CHKERRQ(ierr);

  return ierr;
}

void ElementSet::AssembleTangentStiffness(Mat&A, Vec &B)
{
  AssembleMatrix(A, B, 2, "getTangentStiffness");
}

void ElementSet::AssembleInternalForce(Vec &B)
{
  Mat A;
  AssembleMatrix(A, B, 1, "getInternalForce");
}

void ElementSet::AssembleMassMatrix(Mat &A, Vec &B)
{
  AssembleMatrix(A, B, 2, "getMassMatrix");
}

void ElementSet::CommitHistory()
{
  for(auto elem : m_elem)
    elem.second->CommitHistory();
}

int ElementSet::ElemGroupCount(const std::string&groupName)
{
  if("All" == groupName){
    int temp = 0;
    for(auto group : m_groups)
     temp += group.second.size();
    return temp;
  }
  else
    return m_groups[groupName].size();
}

std::vector<int> ElementSet::IterElementGroup(const std::string&elemGroupName)
{
  std::vector<int> elem;
  if("All" == elemGroupName){
    for(auto iElem : m_groups)
      elem.insert(elem.end(), iElem.second.begin(), iElem.second.end());
  }
  else{
    elem.insert(elem.end(), m_groups[elemGroupName].begin(), m_groups[elemGroupName].end());
  }
  return elem;
}
