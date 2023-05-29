#ifndef NODESET_H
#define NODESET_H

#include <map>

#include <util/Tools.h>

class NodeSet
{
public:
  /**
   * @Brief: Construct a new Node Set object
   * 
   */
  NodeSet();

  /**
   * @Brief: Destroy the Node Set object
   * 
   */
  ~NodeSet();

public:
  /**
   * @Brief: Get the Number Of all Nodes
   * 
   * @return int 
   */
  inline int GetNumOfNodes(){
    return m_nodeCoords.size();
  }
  
  /**
   * @Brief:  Read Node Data form File
   * 
   * @param fileName 
   */
  void ReadFromFile(const std::string&fileName);

  /**
   * @Brief: Get Coordinates of Single Node
   * 
   * @param nodeId                 [in]  Node Index
   * @return std::vector<double>   [out] Node Coordinates
   */
  std::vector<double> GetNodeCoords(const int&nodeId);

  /**
   * @Brief:   Get Coordinates of Multiple Nodes 
   * 
   * @param nodeIds                            [in]  Node Indexes
   * @return std::vector<std::vector<double>>  [out] Node Coordinates
   */
  std::vector<std::vector<double>> GetNodeCoords(const std::vector<int> &nodeIds);

  void UpdateNodeCoords(Vec &dDisp, const int numOfDof);

public:
  std::map<int, std::vector<double>> m_nodeCoords;   // All Node Coordinates
};

#endif // NODESET_H