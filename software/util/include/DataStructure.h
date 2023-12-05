
#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include <petscksp.h>
#include "../../fem/include/NodeSet.h"
#include "../../fem/include/ElementSet.h"
#include "../../fem/include/DofSpace.h"
#include "../../nlohmann/json.hpp"


class GlobalData
{
public:
  static GlobalData* GetInstance();

  static void DestoryInstance();

  void SetFEMData(const nlohmann::json& props, std::shared_ptr<NodeSet> nodes,
                  std::shared_ptr<ElementSet> elements, std::shared_ptr<DofSpace> dofs);

  void ReadFromFile(const std::string&fileName);

  inline void ResetNodalOutput(){
    m_outputName.clear();
  }

  MatrixXd GetData(const std::string&outputName);

  VectorXd GetData(const std::string&outputName, const int nodeID);

  void PrintNodes();

public:
  int m_cycle = 0;
  int m_iiter = 0;
  bool m_active = true;
  double m_time = 0.;
  double m_lam = 0.;
  std::shared_ptr<NodeSet> m_nodes;
  std::shared_ptr<ElementSet> m_elements;
  std::shared_ptr<DofSpace> m_dofs;
  nlohmann::json m_props;
  std::string m_prefix;
  
  std::vector<std::string> m_outputName;
  std::unordered_map<std::string, MatrixXd> m_outputData;
  Vec m_state;
  Vec m_Dstate;
  Vec m_fint;
  Vec m_fhat;
  Vec m_velo;
  Vec m_acce;
private:
  GlobalData();
  ~GlobalData();

  PetscErrorCode CreateVecSpace();

  PetscErrorCode DestroyVecSpace();

  /**
   * @Brief: Read External Nodal Force Information
   * 
   * @param fileName 
   */
  void ReadExternalForce(const std::string&fileName);

  /**
   * @Brief: Read External Nodal Velocity Information
   * 
   * @param fileName 
   */
  void ReadInitialVelocity(const std::string&fileName);

  /**
   * @Brief: Conmon Method to Read External Nodal Force and Velocity
   * 
   * @param data 
   * @param fileName 
   * @param key 
   */
  void ReadData(Vec &data, const std::string&fileName, const std::string &key);


  /**
   * @Brief: Read Edge Loads Information
   * 
   */
  void ReadEdgeLoadsData(const std::string &fileName);

private:
  static GlobalData *m_globalData;
};

#endif // DATASTRUCTURE_H