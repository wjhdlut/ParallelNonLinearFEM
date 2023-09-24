#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H

#include <nlohmann/json.hpp>
#include <materials/BaseMaterial.h>
#include <util/Kinematics.h>


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
   * @param iSam 
   * @return std::vector<double> 
   */
  VectorXd GetStress(const std::shared_ptr<Kinematics>&kin,
                     const VectorXd &increDisp = VectorXd::Zero(0),
                     const MatrixXd &dphi = MatrixXd::Zero(0, 0), int iSam = -1);
                     
  /**
   * @Brief:  Get the Tangent Modulue Matrix
   * 
   * @return std::vector<std::vector<double>> 
   */
  inline MatrixXd GetTangMatrix(){
    return m_mat->GetTangMatrix();
  }

  inline MatrixXd GetTangMatrix(const MatrixXd &F){
    return m_mat->GetTangentMatrix(F);
  }

  /**
   * @Brief:         
   * 
   */
  inline void CommitHistory(){
    m_mat->CommitHistory();
  }

  inline double GetMaterialRho(){
    return m_mat->ReturnRho();
  }

  inline double GetMaterialPara(const std::string &name){
    return m_mat->SetMaterialParamter(name);
  }

private:
  int iIter = -1;                        // 
  std::shared_ptr<BaseMaterial> m_mat;   // Material Pointer
};

#endif // MATERIALMANAGER_H

