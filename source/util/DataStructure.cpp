#include <fstream>
#include <iostream>
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
        if(line.npos != line.find("</ExternalForces>")) return;

        std::string temp = Tools::StringStrip(line);
        std::vector<std::string> a = Tools::StringSplit(temp, ";");

        std::vector<std::string> b = Tools::StringSplit(a[0], "=");
        std::vector<std::string> c = Tools::StringSplit(b[0], "[");

        std::string dofType = c[0];
        int nodeId = std::stoi(Tools::StringSplit(c[1], "]")[0]);
        std::cout << m_dofs->GetForType(nodeId, dofType) << std::endl;

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