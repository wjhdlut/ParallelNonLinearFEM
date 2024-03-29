
#include <fstream>
#include <regex>
#include <iostream>
#include <fem/NodeSet.h>
#include <util/DataStructure.h>


NodeSet::NodeSet()
{
  m_nodeCoords = {};
}

NodeSet::~NodeSet()
{}

void NodeSet::ReadFromFile(const std::string&fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  if(!fin.is_open()){
    std::cout << fileName << " open failed!!" << std::endl;
    exit(-1);
  }
  std::string line = "";

  std::regex pattern("[\\s]{2,}");
  std::regex intPattern();
  std::smatch result;
  while(true)
  {
    getline(fin , line);
    if(line.npos != line.find("\r"))
      line.erase(line.find("\r"));
    if(line.npos != line.find("<Nodes>"))
    {
      while(true)
      {
        getline(fin, line);
        if(line.npos != line.find("\r"))
          line.erase(line.find("\r"));

        if(line.npos != line.find("</Nodes>")) return;

        if(line.size() == 0) continue; 
        
        line =  std::regex_replace(line, pattern, " ");
        std::vector<std::string> strData = Tools::StringSplit(line, ";");
        for(auto iter = strData.begin(); iter != strData.end() - 1; iter++)
        {
          std::string tempStr = Tools::StringStrip(*iter);
          std::vector<std::string> b = Tools::StringSplit(tempStr, " ");
          if("//" == b[0].substr(0, 2) || "#" == b[0].substr(0, 1)) break;

          int nodeID = std::stoi(b[0]);

          std::vector<double> coords;
          for(auto iterb = b.begin() + 1; iterb != b.end(); iterb++)
            coords.emplace_back(std::stod(*iterb));
          
          if(m_nodeCoords.count(nodeID) == 0){
            m_nodeCoords[nodeID] = VectorXd::Zero(coords.size());
            for(unsigned i = 0; i < coords.size(); i++)
              m_nodeCoords[nodeID](i) = coords[i];
          }
          
        }
      }
    }
  }
  fin.close();
}

VectorXd NodeSet::GetNodeCoords(const int&nodeId)
{
  if(0 == m_nodeCoords.count(nodeId)){
    throw "Node ID " + std::to_string(nodeId) + " does not exist";
  }
  return m_nodeCoords.at(nodeId);
}

MatrixXd NodeSet::GetNodeCoords(const std::vector<int> &nodeIds)
{
  MatrixXd coords = MatrixXd::Zero(nodeIds.size(), m_nodeCoords[nodeIds[0]].size());
  for(int i = 0; i < nodeIds.size(); i++)
    coords.row(i) = GetNodeCoords(nodeIds[i]);
  return coords;
}

void NodeSet::UpdateNodeCoords(Vec&dDisp, const int numOfDof)
{
  std::vector<int> indexDof;
  VectorXd iNodeDDisp = VectorXd::Zero(numOfDof);
  for(auto iNode = m_nodeCoords.begin(); iNode != m_nodeCoords.end(); iNode++)
  {
    indexDof = GlobalData::GetInstance()->m_dofs->GetForType(iNode->first);
    VecGetValues(dDisp, numOfDof, &indexDof[0], &iNodeDDisp(0));
    iNode->second += iNodeDDisp;
  }
}