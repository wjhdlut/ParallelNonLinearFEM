#include <fstream>
#include <regex>
#include <iostream>
#include <iomanip>
#include "../include/DataStructure.h"

GlobalData *GlobalData::m_globalData = nullptr;

GlobalData::GlobalData()
{
}

GlobalData::~GlobalData()
{
  DestroyVecSpace();
}

GlobalData* GlobalData::GetInstance()
{
  if(nullptr == m_globalData){
    m_globalData = new GlobalData();
  }

  return m_globalData;
}

void GlobalData::DestoryInstance()
{
  if(nullptr != m_globalData) {
    delete m_globalData;
  }
}

void GlobalData::SetFEMData(const nlohmann::json& props, std::shared_ptr<NodeSet> nodes,
                            std::shared_ptr<ElementSet> elements, std::shared_ptr<DofSpace> dofs)
{
  m_props = props;
  m_nodes = nodes;
  m_elements = elements;
  m_dofs = dofs;
  CreateVecSpace();
}

void GlobalData::ReadFromFile(const std::string&fileName)
{
  ReadExternalForce(fileName);

  ReadInitialVelocity(fileName);

  ReadEdgeLoadsData(fileName);
}

PetscErrorCode GlobalData::CreateVecSpace()
{
  PetscErrorCode ierr;
  int numOfDofs = m_dofs->m_dofs.size() * m_dofs->m_dofs.at(0).size();
  ierr = VecCreate(PETSC_COMM_WORLD, &m_state); CHKERRQ(ierr);
  ierr = VecSetSizes(m_state, PETSC_DECIDE, numOfDofs); CHKERRQ(ierr);
  ierr = VecSetFromOptions(m_state); CHKERRQ(ierr);
  ierr = VecSet(m_state, 0.0); CHKERRQ(ierr);

  ierr = VecDuplicate(m_state, &m_Dstate); CHKERRQ(ierr);
  // ierr = VecCopy(m_state, m_Dstate); CHKERRQ(ierr);
  
  ierr = VecDuplicate(m_state, &m_fint); CHKERRQ(ierr);
  // ierr = VecCopy(m_state, m_fint); CHKERRQ(ierr);
  
  ierr = VecDuplicate(m_state, &m_fhat); CHKERRQ(ierr);
  // ierr = VecCopy(m_state, m_fhat); CHKERRQ(ierr);
  
  ierr = VecDuplicate(m_state, &m_velo); CHKERRQ(ierr);
  // ierr = VecCopy(m_state, m_velo); CHKERRQ(ierr);
  
  ierr = VecDuplicate(m_state, &m_acce); CHKERRQ(ierr);
  // ierr = VecCopy(m_state, m_acce); CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode GlobalData::DestroyVecSpace()
{
  PetscErrorCode ierr;
  ierr = VecDestroy(&m_state); CHKERRQ(ierr);
  ierr = VecDestroy(&m_Dstate); CHKERRQ(ierr);
  ierr = VecDestroy(&m_fint); CHKERRQ(ierr);
  ierr = VecDestroy(&m_fhat); CHKERRQ(ierr);
  ierr = VecDestroy(&m_velo); CHKERRQ(ierr);
  ierr = VecDestroy(&m_acce); CHKERRQ(ierr);
  return ierr;
}

MatrixXd GlobalData::GetData(const std::string&outputName)
{
  MatrixXd data    = m_outputData[outputName];
  MatrixXd weights = m_outputData[outputName + "Weight"];

  MatrixXd outputData = MatrixXd::Zero(m_nodes->m_nodeCoords.size(), data.cols());
  VectorXd iRow;
  
  int count = 0;
  for(auto node : m_nodes->m_nodeCoords)
  {
    int index = m_dofs->GetIndex(node.first);
    if(weights(index, 0) != 0)
      iRow = 1./weights(index, 0) * data.row(index);
    else
      iRow = data.row(index);

    outputData.row(count) = iRow;

    count += 1;
  }

  return outputData;
}

VectorXd GlobalData::GetData(const std::string&outputName, const int nodeID)
{
  MatrixXd data = m_outputData[outputName];
  MatrixXd weights = m_outputData[outputName + "Weight"];
  
  int index = m_dofs->GetIndex(nodeID);
  int weightValue = weights(index, 0);
  
  VectorXd tempVec(data.row(index).size());
  // for(auto iData : data[index])
  //   (0 != weightValue) ? tempVec.emplace_back(iData/weightValue) : tempVec.emplace_back(iData);
    
  for(int i = 0; i < data.row(index).size(); i++)
    tempVec(i) = (0 != weightValue) ? data(index, i) / weightValue : data(index, i);
  return tempVec;
}

void GlobalData::PrintNodes()
{
  std::cout << "   Node | ";

  for(auto dofType : m_dofs->GetDofType())
    std::cout << "    " << dofType << "         ";

  int intTemp;
  VecGetSize(m_fint, &intTemp);
  if(0 != intTemp)
  {
    for(auto dofType : m_dofs->GetDofType())
      std::cout << "  fint-" << dofType << "      ";
  }

  for(auto name : m_outputName)
    std::cout << "         " << name;
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------"
            << "-------------------------------------------------------\n";
  
  std::vector<double> doubleTemp(2*m_dofs->GetDofType().size());
  std::vector<int> index;
  for(auto iNode : m_nodes->m_nodeCoords)
  {
    std::cout << std::setw(6) << std::right << iNode.first << "  |";
    
    // Output the Displacement
    index.clear();
    for(auto dofType : m_dofs->GetDofType())
      index.emplace_back(m_dofs->GetForType(iNode.first, dofType));
    
    VecGetValues(m_state, index.size(), &index[0], &doubleTemp[0]);
    VecGetValues(m_fint, index.size(), &index[0], &doubleTemp[index.size()]);

    // Output the Fint Values
    std::cout.precision(3);
    for(auto iTempValue : doubleTemp)
      std::cout << std::setw(13) << std::setiosflags(std::ios::scientific) << iTempValue;

    for(auto name : m_outputName)
    {
      VectorXd data = GetData(name, iNode.first);
      for(int i = 0; i < data.size(); i++)
        std::cout << std::setw(13) << std::setiosflags(std::ios::scientific) << data(i);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void GlobalData::ReadData(Vec &data, const std::string&fileName, const std::string &key)
{
  std::ifstream fin(fileName, std::ios::in);
  std::string line = "";

  while(true){
    getline(fin, line);
    if(line.npos != line.find("\r")) line.erase(line.find("\r"));

    if(line.npos != line.find(('<' + key + '>')))
    {
      while(true){
        getline(fin, line);
        if(line.npos != line.find("\r")) line.erase(line.find("\r"));
        line.erase(0, line.find_first_not_of(" "));
        if(0 == line.size()) continue;
        if(line.npos != line.find(("</" + key + '>'))){
          VecAssemblyBegin(data);
          VecAssemblyEnd(data);
          fin.close();
          return;
        }

        if(0 == line.size()) continue;

        std::string temp = Tools::StringStrip(line);
        std::vector<std::string> a = Tools::StringSplit(temp, ";");

        std::vector<std::string> b = Tools::StringSplit(a[0], "=");
        std::vector<std::string> c = Tools::StringSplit(b[0], "[");

        std::string dofType = c[0];
        int nodeId = std::stoi(Tools::StringSplit(c[1], "]")[0]);

        VecSetValue(data, m_dofs->GetForType(nodeId, dofType), std::stod(b[1]), ADD_VALUES);
      }
    }
    if(fin.eof()) break;
  }
  fin.close();
}

void GlobalData::ReadExternalForce(const std::string&fileName)
{
  ReadData(m_fhat, fileName, "ExternalForces");
}

void GlobalData::ReadInitialVelocity(const std::string&fileName)
{
  ReadData(m_velo, fileName, "InitialVelocity");
}

void GlobalData::ReadEdgeLoadsData(const std::string &fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  std::string line = "";
  std::regex pattern("[\\s]{2,}");

  while(true){
    getline(fin, line);

    if(line.npos != line.find("\r")) line.erase(line.find("\r"));

    if(line.npos != line.find(("<EdgeLoads>")))
    {
      while(true){
        getline(fin, line);
        if(line.npos != line.find("\r")) line.erase(line.find("\r"));

        line.erase(0, line.find_first_not_of(" "));

        if(0 == line.size()) continue;

        if(line.npos != line.find("</EdgeLoads>")){
          VecAssemblyBegin(m_fhat);
          VecAssemblyEnd(m_fhat);
          fin.close();
          return;
        }

        line =  std::regex_replace(line, pattern, " ");
        std::vector<std::string> a = Tools::StringSplit(line, ";");
        std::vector<std::string> b = Tools::StringSplit(a[0], " ");
        
        int elemID = std::stoi(b[0]);
        int numOfNode = std::stoi(b[1]);
        
        // read edge load information
        std::shared_ptr<Element> elemPtr = m_elements->GetElementPtr()[elemID];
        std::unordered_map<int, std::vector<double>> nodeForcePres;
        std::unordered_map<int, VectorXd> nodeCoords;
        
        for(int iNode = 0; iNode < numOfNode; iNode++){
          int nodeId = std::stoi(b[2+iNode]);
          nodeForcePres.insert(std::pair<int, std::vector<double>>(nodeId, {}));
          nodeCoords.insert(std::pair<int, VectorXd>(nodeId, m_nodes->GetNodeCoords(nodeId)));
          for(int iDof = 0; iDof < m_nodes->GetNodeCoords(nodeId).size(); iDof++)
            nodeForcePres[nodeId].emplace_back(std::stod(b[2 + numOfNode + iDof*numOfNode + iNode]));
        }
        
        // Compute Equivalent Node Force
        VectorXd elemForce = elemPtr->CompEquivalentNodeForce(nodeCoords, nodeForcePres);
        // Assemble into total force vector
        std::vector<int> elemDofs = m_dofs->Get(elemPtr->GetNodes());
        VecSetValues(m_fhat, elemDofs.size(), &elemDofs[0],&elemForce(0), ADD_VALUES);

      }
    }
    if(fin.eof()) break;
  }
  fin.close();
}
                                