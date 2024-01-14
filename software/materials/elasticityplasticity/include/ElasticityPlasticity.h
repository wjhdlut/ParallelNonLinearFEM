/**
 * @File Name:     ElasticityPlasticity.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-26
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef ELASTICITYPLASTICITY_H
#define ELASTICITYPLASTICITY_H

#include "YieldRule.h"
#include "../../include/BaseMaterial.h"
#include "../../linearelasticity/include/LinearElasticity.h"
#include "../../../util/include/Math.h"
#include "../../../util/include/ObjectFactory.h"

class ElasticityPlasticity : public BaseMaterial
{
public:
  /**
   * @Brief: Construct a new Elasticity Plasticity object
   * 
   * @param matProps 
   */
  ElasticityPlasticity(const nlohmann::json &matProps);

  /**
   * @Brief: Destroy the Elasticity Plasticity object
   * 
   */
  ~ElasticityPlasticity();

  /**
   * @Brief: Get the Stress Vector
   * 
   * @param kin 
   * @param increDisp 
   * @param dphi 
   * @return VectorXd 
   */
  virtual VectorXd GetStress(const std::shared_ptr<Kinematics> &kin,
                             const VectorXd &stress = VectorXd::Zero(0)) override;

protected:
  virtual void ComputeDMatrix() override;

  /**
   * @Brief: Rotation Stress to Get Objective Stress Rate
   * 
   * @param stress 
   * @param increDisp 
   * @param dphi 
   */
  void StressRotation(VectorXd &stress,
                      const VectorXd &increDisp,
                      const MatrixXd &dphi);

  /**
   * @Brief: Initialize Some Variable;
   * 
   */
  void Initialize();

private:
  /**
   * @Brief: Read Yield Stress and Strain Data
   * 
   */
  void ReadYieldStressStrainData();
  
  /**
   * @Brief: Get the Yield Stress for the Given Accumulated Plastic Strain
   *         Based on the Stress Strain Data Pair
   * 
   * @param plasticStrain 
   * @return double 
   */
  double GetYieldStress(const double accumPlasticStrain);

   /**
   * @Brief: Get the Yield Stress for the Given Accumulated Plastic Strain
   *         Based on the Stress Strain Data Pair
   * 
   * @param plasticStrain 
   * @return double 
   */
  double GetHardModuli(const double accumPlasticStrain);

protected:
  bool m_yieldFlag = false;                                       // Plastic Yielding Flag
  double m_K;                                                     // Bulk Moduli
  double m_G;                                                     // Shear Moduli
  double m_plaMod;                                                // Hard moduli
  double m_waveSpeed;                                             // Wave Speed for Explicit Dynamic Problem 
  double m_accumPlasticStrain;                                    // Accumulated Plastic Strain
  double m_plasticMulter;
  std::string m_rateType;                                         // Type of Objective Rate
  VectorXd m_oneVec;                                              // Common Vector
  VectorXd m_plasticStrain;                                       // Plastic Strain
  MatrixXd m_devMat;                                              // Common Matrix
  MatrixXd m_oneMat;                                              // Common Matrix
  std::vector<std::pair<double, double>> m_yieldStress;           // Harding Curve
  std::shared_ptr<LinearElasticity> m_lineMat;                    // Instance of Linear Elastic Matrial
  std::shared_ptr<YieldRule> m_yieldRule;                         // Yield Rule

private:
  double m_tol = 1.0e-6;
};
ReflectRegister(ElasticityPlasticity, const nlohmann::json &)

#endif // ELASTICITYPLASTICITY_H