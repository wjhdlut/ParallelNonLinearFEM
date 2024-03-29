
#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include <fem/NodeSet.h>
#include <fem/ElementSet.h>
#include <fem/DofSpace.h>
#include <nlohmann/json.hpp>
#include <petscksp.h>


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

  void ReadExternalForce(const std::string&fileName);

  void ReadInitialVelocity(const std::string&fileName);

  void ReadData(Vec &data, const std::string&fileName, const std::string &key);

private:
  static GlobalData *m_globalData;
};

#endif // DATASTRUCTURE_H