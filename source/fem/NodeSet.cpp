
#include <fstream>
#include <regex>
#include <iostream>
#include <fem/NodeSet.h>


NodeSet::NodeSet()
{
  m_nodeCoords = {};
}

NodeSet::~NodeSet()
{}

void NodeSet::ReadFromFile(const std::string&fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  std::string line = "";

  std::regex pattern("[\\s]{2,}");
  std::regex intPattern();
  std::smatch result;
  while(true)
  {
    getline(fin , line), line.erase(line.find("\r"));
    if(line.npos != line.find("<Nodes>"))
    {
      while(true)
      {
        getline(fin, line), line.erase(line.find("\r"));

        if(line.npos != line.find("</Nodes>")) return;
        
        line =  std::regex_replace(line, pattern, " ");
        std::vector<std::string> strData = Tools::StringSplit(line, ";");
        for(auto iter = strData.begin(); iter != strData.end() - 1; iter++)
        {
          std::string tempStr = Tools::StringStrip(*iter);
          std::vector<std::string> b = Tools::StringSplit(tempStr, " ");
          if("//" == b[0].substr(0, 2) || "#" == b[0].substr(0, 1)) break;

          int nodeID = std::stoi(b[0]);
          // m_nodeIndex.emplace_back(nodeID);
          if(m_nodeCoords.count(nodeID) == 0)
            m_nodeCoords.insert(std::pair<int, std::vector<double>>(nodeID, {}));
          for(auto iterb = b.begin() + 1; iterb != b.end(); iterb++)
            m_nodeCoords[nodeID].emplace_back(std::stod(*iterb));
        }
      }
    }
  }
  fin.close();
}

std::vector<double> NodeSet::GetNodeCoords(const int&nodeId)
{
  if(0 == m_nodeCoords.count(nodeId)){
    throw "Node ID " + std::to_string(nodeId) + " does not exist";
  }
  return m_nodeCoords.at(nodeId);
}

std::vector<std::vector<double>> NodeSet::GetNodeCoords(const std::vector<int> &nodeIds)
{
  std::vector<std::vector<double>> coords;
  for(auto iNode : nodeIds)
    coords.emplace_back(GetNodeCoords(iNode));
  return coords;
}