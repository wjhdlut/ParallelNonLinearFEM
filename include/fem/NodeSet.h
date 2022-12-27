#ifndef NODESET_H
#define NODESET_H

#include <map>

#include <util/Tools.h>

class NodeSet
{
public:
  NodeSet();
  ~NodeSet();

public:
  /**
   * @Brief:         Get the Number Of Nodes
   * 
   * @return int 
   */
  inline int GetNumOfNodes(){
    return m_nodeCoords.size();
  }

  void ReadFromFile(const std::string&fileName);

  std::vector<double> GetNodeCoords(const int&nodeId);

  std::vector<std::vector<double>> GetNodeCoords(const std::vector<int> &nodeIds);

public:
  std::map<int, std::vector<double>> m_nodeCoords;
};

#endif // NODESET_H