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


void DofSpace::Solve(Mat&K, Vec&df, Vec&da)
{

}

void DofSpace::Solve(double K, Vec&df, Vec&da)
{
  if (abs(K) > 1e-10)
  {
    Vec temp;
    VecDuplicate(da, &temp);
    VecSet(temp, 0.);
    VecWAXPY(da, 1. / K, df, temp);
    VecDestroy(&temp);
  }
}

Matrix DofSpace::GetConstraintsMatrix()
{
  int n_constrianed = m_constrained.size();
  int n = m_dofs.size();

  Matrix C = Math::MatrixZeros(n, n-n_constrianed);
  
  int j = 0;
  for(int i = 0; i < n ; i++)
  {
    if(0 != m_constrained.count(i)) continue;
    C[i][j] = 1.;
    j += 1;
  }

  return C;
}