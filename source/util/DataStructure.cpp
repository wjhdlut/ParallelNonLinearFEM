#include <fstream>
#include <iostream>
#include <iomanip>
#include <util/DataStructure.h>

GlobalData *GlobalData::m_globalData = nullptr;

GlobalData::GlobalData()
{
}

GlobalData::~GlobalData()
{
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
  std::ifstream fin(fileName, std::ios::in);

  std::string line = "";

  while(true){
    getline(fin, line), line.erase(line.find("\r"));

    if(line.npos != line.find("<ExternalForces>"))
    {
      while(true){
        getline(fin, line), line.erase(line.find("\r"));
        line.erase(0, line.find_first_not_of(""));
        if(0 == line.size()) continue;
        if(line.npos != line.find("</ExternalForces>")) return;

        std::string temp = Tools::StringStrip(line);
        std::vector<std::string> a = Tools::StringSplit(temp, ";");

        std::vector<std::string> b = Tools::StringSplit(a[0], "=");
        std::vector<std::string> c = Tools::StringSplit(b[0], "[");

        std::string dofType = c[0];
        int nodeId = std::stoi(Tools::StringSplit(c[1], "]")[0]);

        VecSetValue(m_fhat, m_dofs->GetForType(nodeId, dofType), std::stod(b[1]), ADD_VALUES);
      }
    }
  }
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
  ierr = VecDuplicate(m_state, &m_fint); CHKERRQ(ierr);
  ierr = VecDuplicate(m_state, &m_fhat); CHKERRQ(ierr);
  ierr = VecDuplicate(m_state, &m_velo); CHKERRQ(ierr);
  ierr = VecDuplicate(m_state, &m_acce); CHKERRQ(ierr);
  return ierr;
}

Matrix GlobalData::GetData(const std::string&outputName)
{
  Matrix data    = m_outputData[outputName];
  Matrix weights = m_outputData[outputName + "Weight"];

  Matrix outputData;
  std::vector<double> iRow;
  
  for(auto node : m_nodes->m_nodeCoords)
  {
    int index = m_dofs->GetIndex(node.first);
    if(weights[index][0] != 0)
      iRow = Math::VecScale(1./weights[index][0], data[index]);
    else
      iRow = data[index];

    outputData.emplace_back(iRow);
  }

  return outputData;
}

std::vector<double> GlobalData::GetData(const std::string&outputName, const int nodeID)
{
  Matrix data = m_outputData[outputName];
  Matrix weights = m_outputData[outputName + "Weight"];
  
  int index = m_dofs->GetIndex(nodeID);
  int weightValue = weights[index][0];
  std::vector<double> tempVec;
  for(auto iData : data[index])
    (0 != weightValue) ? tempVec.emplace_back(iData/weightValue) : tempVec.emplace_back(iData);
    
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
      auto data = GetData(name, iNode.first);
      for(auto iData : data)
        std::cout << std::setw(13) << std::setiosflags(std::ios::scientific) << iData;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}