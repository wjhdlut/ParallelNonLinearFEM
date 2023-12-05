#ifndef NODESET_H
#define NODESET_H

#include <map>

#include "../../util/include/Tools.h"
#include <eigen3/Eigen/Dense>

using namespace Eigen;

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
  VectorXd GetNodeCoords(const int&nodeId);

  /**
   * @Brief:   Get Coordinates of Multiple Nodes 
   * 
   * @param nodeIds                            [in]  Node Indexes
   * @return std::vector<std::vector<double>>  [out] Node Coordinates
   */
  MatrixXd GetNodeCoords(const std::vector<int> &nodeIds);

  /**
   * @Brief: Update Node Coordinates
   * 
   * @param dDisp 
   * @param numOfDof 
   */
  void UpdateNodeCoords(Vec &dDisp, const int numOfDof);

public:
  std::map<int, VectorXd> m_nodeCoords;   // All Node Coordinates

private:
  std::string m_coordSysType = "Rectangular";  // Type of Coordinate system

private:
  VectorXd TransCoordCylToRec(const std::vector<double> &coord);
  const double Pi = 3.1415926;
};

#endif // NODESET_H