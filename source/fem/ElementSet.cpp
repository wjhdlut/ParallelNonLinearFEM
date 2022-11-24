#include <fstream>
#include <regex>
#include <fem/ElementSet.h>
#include <elements/Element.h>
#include <util/ObjectFactory.h>

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
  std::string line = "";
  std::regex pattern("[\\s]{2,}");
  while(true)
  {
    getline(fin, line); line.erase(line.find("\r"));

    if(line.npos != line.find("<Elements>"))
    {
      while(true)
      {
        getline(fin, line); line.erase(line.find("\r"));

        if(line.npos != line.find("</Elements>")) return;
        
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
          
          Add(elemID, b[1], elementNodes);
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

  for(auto nodeID : elem->GetNodes())
    m_nodes->GetNodeCoords(nodeID);

  m_elem[elemId] = elem;
  if (0 == m_groups.count(modelName)) m_groups[modelName] = {};
  m_groups[modelName].emplace_back(elemId);
}