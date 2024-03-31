/**
 * @File Name:     MaterialManager.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-26
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H

#include "BaseMaterial.h"
#include "../../nlohmann/json.hpp"
#include "../../util/include/Kinematics.h"


/**
 * @Brief: class of material manager to manage materials of model
 * 
 */
class MaterialManager
{
public:
  /**
  * @Brief: Construct a new Material Manager object
  * 
  * @param matProps      [in] The Properties of Whole Models
  */
  MaterialManager(const nlohmann::json &matProps);

  /**
   * @Brief: Destroy the Material Manager object
   * 
   */
  ~MaterialManager();

  /**
   * @Brief:  Reset Member Varibale <iIter> data
   * 
   */
  void Reset();
  
  /**
   * @Brief:  Get the Stress Data
   * 
   * @param kin
   * @param stress     [in] Stress at the Last Load Step
   * @param iSam 
   * @return std::vector<double> 
   */
  VectorXd GetStress(const std::shared_ptr<Kinematics>&kin,
                     int iSam = -1);
                     
  /**
   * @Brief:  Get the Tangent Modulue Matrix
   * 
   * @return std::vector<std::vector<double>> 
   */
  inline MatrixXd GetTangMatrix(){
    return m_mat->GetTangMatrix();
  }
  
  /**
   * @Brief:  Get the Tang Matrix in Initial Configuration
   * 
   * @param F 
   * @return MatrixXd 
   */
  inline MatrixXd GetTangMatrix(const MatrixXd &F){
    return m_mat->GetTangentMatrix(F);
  }

  /**
   * @Brief: Commit History 
   * 
   */
  inline void CommitHistory(){
    m_mat->CommitHistory();
  }

  /**
   * @Brief: Get the Material Density
   * 
   * @return double 
   */
  inline double GetMaterialRho(){
    return m_mat->ReturnRho();
  }
  
  /**
   * @Brief: Get the Material Parameter
   * 
   * @param name 
   * @return double 
   */
  inline double GetMaterialPara(const std::string &name){
    return m_mat->SetMaterialParamter(name);
  }

  // /**
  //  * @Brief: Get the Element State Variable
  //  * 
  //  * @param stateVariable 
  //  */
  // inline void GetElemStateVariable(const double &stateVariable){
  //   m_mat->GetElemStateVariable(stateVariable);
  // }

  // inline void SetElemStateVariable(double &stateVariable){
  //   m_mat->SetElemStateVariable(stateVariable);
  // }

private:
  int iIter = -1;                        // 
  std::shared_ptr<BaseMaterial> m_mat;   // Material Pointer
};

#endif // MATERIALMANAGER_H

