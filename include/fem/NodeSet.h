#ifndef NODESET_H
#define NODESET_H

#include <unordered_map>

#include <util/Tools.h>

class NodeSet
{
public:
  NodeSet();
  ~NodeSet();

public:
  void ReadFromFile(const std::string&fileName);

  std::vector<double> GetNodeCoords(const int&nodeId);

  std::vector<std::vector<double>> GetNodeCoords(const std::vector<int> &nodeIds);

private:
  std::unordered_map<int, std::vector<double>> m_nodeCoords;
};

#endif // NODESET_H