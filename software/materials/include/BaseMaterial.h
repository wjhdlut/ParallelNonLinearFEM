/**
 * @File Name:     BaseMaterial.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-26
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef BASEMATERIAL_H
#define BASEMATERIAL_H

#include <vector>
#include <unordered_map>

#include "../../nlohmann/json.hpp"
#include "../../util/include/Kinematics.h"

class BaseMaterial
{
public:
  /**
   * @Brief: Construct a new Base Material object
   * 
   * @param props 
   */
  BaseMaterial(const nlohmann::json &props);

  /**
   * @Brief: Destroy the Base Material object
   * 
   */
  virtual ~BaseMaterial();

  /**
   * @Brief: Get the Stress Vector
   * 
   * @param kin 
   * @param stress
   * @return VectorXd 
   */
  virtual VectorXd GetStress(const std::shared_ptr<Kinematics>&kin,
                             const VectorXd &stress = VectorXd::Zero(0));

  /**
   * @Brief: Set the Iter object
   * 
   * @param iIter 
   */
  void SetIter(int iIter);

  /**
   * @Brief: Commit History Information
   * 
   */
  void CommitHistory();

  /**
   * @Brief: Set the Material Paramter
   * 
   * @param name 
   * @return double 
   */
  double SetMaterialParamter(const std::string &name);

  /**
   * @Brief: Get the Tangent Matrix
   * 
   * @return MatrixXd 
   */
  inline MatrixXd GetTangMatrix(){
    if(0 == m_D.rows()) ComputeDMatrix();
    if(m_updateDMatrix) ComputeDMatrix();
    return m_D;
  }
  
  /**
   * @Brief: Pull Back the Tangent Matrix to Initial Configuration
   * 
   * @param F 
   * @return MatrixXd 
   */
  inline MatrixXd GetTangentMatrix(const MatrixXd &F){
    if(0 == m_D.rows()) ComputeDMatrix();
    TangentDMatrixToInitial(F);
    return m_D;
  }
  
  /**
   * @Brief: Return Material Density
   * 
   * @return double 
   */
  inline double ReturnRho(){
    return m_rho;
  }

protected:
  /**
   * @Brief: Set the History Parameter
   * 
   * @param name 
   * @param value 
   * @return int 
   */
  int SetHistoryParameter(const std::string &name, double value);

  /**
   * @Brief: Get the History Parameter
   * 
   * @param name 
   * @return double 
   */
  double GetHistoryParameter(const std::string&name);

  /**
   * @Brief: Compute Tangent Matrix
   * 
   */
  virtual void ComputeDMatrix() = 0;

  /**
   * @Brief: Push Forward Tangent Matrix to Current Configuration
   * 
   * @param F 
   */
  void TangentDMatrixToCurrent(const MatrixXd &F);

  /**
   * @Brief: Pull Back the Tangent Matrix to Initial Configuration
   * 
   * @param F 
   */
  void TangentDMatrixToInitial(const MatrixXd &F);

private:
  /**
   * @Brief: Get the Transform Matrix
   * 
   * @param T 
   * @param F 
   */
  void GetTransMatrix(MatrixXd &T, const MatrixXd &F);

  /**
   * @Brief: Initialize Some Basic Variables
   * 
   */
  void Initialize();

protected:
  bool   m_planeStrainFlag = false;
  bool   m_planeStressFlag = false;
  bool   m_updateDMatrix   = false;
  int    m_iIter  = -1;
  double m_rho    = 0.;
  double m_E      = 0.;
  double m_nu     = 0.;
  std::vector<std::unordered_map<std::string, double>> m_current;
  std::vector<std::unordered_map<std::string, double>> m_history;
  std::unordered_map<std::string, double> m_initHistory;
  nlohmann::json m_props;
  MatrixXd m_D;
};

#endif // BASEMATERIAL_H