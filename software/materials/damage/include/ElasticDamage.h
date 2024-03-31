/**
 * @File Name:     ElasticDamage.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2024-01-21
 * 
 * @Copyright Copyright (c) 2024 JianHuaWang
 * 
 */

#ifndef ELASTICDAMAGE_H
#define ELASTICDAMAGE_H

#include "EquivantStrainBase.h"
#include "../../include/BaseMaterial.h"
#include "../../../util/include/ObjectFactory.h"
#include "../../linearelasticity/include/LinearElasticity.h"

class ElasticDamage : public BaseMaterial
{
public:
  /**
   * @Brief: Default Constructor of ElasticDamage Class
   * 
   */
  ElasticDamage() = default;
  
  /**
   * @Brief: Constructor of ElasticDamage Class
   * 
   * @param props 
   */
  ElasticDamage(const nlohmann::json &props);

  /**
   * @Brief: Destructor of ElasticDamage Class
   * 
   */
  ~ElasticDamage();

public:
  /**
   * @Brief:  Get the Stress
   * 
   * @param kin 
   * @param stress 
   * @return VectorXd 
   */
  virtual VectorXd GetStress(const std::shared_ptr<Kinematics>&kin);

  virtual void ComputeDMatrix();

private:
  /**
   * @Brief: Initialize Some Basic Variables
   * 
   */
  void Initialize();
  
  /**
   * @Brief: Compute Equivance Strain
   * 
   * @param kin 
   */
  void CompEquivanceStrain(double &eps, 
                           VectorXd &dEpsdStrain,
                           const std::shared_ptr<Kinematics>&kin);

  /**
   * @Brief: Compute the Damage Variable
   * 
   * @param kappa 
   */
  void GetDamage(const double &kappa);
  
  /**
   * @Brief: Compute the Tangent Matrix
   * 
   */
  void CompTangMatrix();

private:
  bool   m_progDam = false;
  double m_sc      = 1./3.;
  double m_kappa   = 0.;
  double m_kappac  = 0.;
  double m_omega   = 0.;
  double m_dOmegadKappa = 0.;
  std::shared_ptr<LinearElasticity> m_linearMat         = nullptr;
  std::shared_ptr<EquivantStrainBase> m_equivanceStrain = nullptr;

  VectorXd m_effStress;
  VectorXd m_dEpsdStrain;
};

ReflectRegister(ElasticDamage, const nlohmann::json &)

#endif // ELASTICDAMAGE_H