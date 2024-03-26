/**
 * @File Name:     ElasticDamage.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2024-03-02
 * 
 * @Copyright Copyright (c) 2024 JianHuaWang
 * 
 */

#include <iostream>

#include "../include/ElasticDamage.h"
#include "../../../util/include/Math.h"
#include "../../../util/include/ObjectFactory.h"

ElasticDamage::ElasticDamage(const nlohmann::json &props) : BaseMaterial(props)
{
  Initialize();
}

ElasticDamage::~ElasticDamage()
{
}

void ElasticDamage::Initialize()
{
  m_linearMat = std::make_shared<LinearElasticity>(m_props);

  std::string damageModel = "VonMise";
  if (m_props.contains("damageModel"))
    damageModel = m_props.at("damageModel");
  m_equivanceStrain = ObjectFactory::CreateObject<EquivantStrainBase>(
                      damageModel + "EquivantStrain", m_props);
  if(nullptr == m_equivanceStrain){
    std::cout << damageModel << " Damage Model Created Failed!!!" << std::endl;
    exit(-1);
  }
  
  // Read Some Basic Parameters
  m_kappa = SetMaterialParamter("kappa0");
  m_kappac = SetMaterialParamter("kappac");

  SetHistoryParameter("kappa", 0.);

  CommitHistory();
}

VectorXd ElasticDamage::GetStress(const std::shared_ptr<Kinematics>&kin,
                                  const VectorXd &stress)
{
  double kappa = GetHistoryParameter("kappa");
  
  // Compute the Effective Strain and its derivative w.r.t. strain array
  double eps = 0;
  // Vector3d m_dEpsdStrain = Vector3d::Zero();
  m_equivanceStrain->CompEquivanceStrain(eps, m_dEpsdStrain, kin);
  // std::cout << "m_dEpsdStrain = \n" << m_dEpsdStrain << std::endl;

  if(eps > kappa){
    m_progDam = true;
    kappa = eps;
    std::cout << "============= Damage =============" << std::endl;
  }
  else{
    m_progDam = false;
  }

  SetHistoryParameter("kappa", kappa);

  // Compute the Damage Variable
  GetDamage(kappa);
  
  // Compute the Effective Stress
  m_effStress = m_linearMat->GetTangMatrix() * kin->strain;
  // std::cout << "m_effStress = \n" << m_effStress << std::endl;
  // std::cout << "m_omega = " << m_omega << std::endl;
  return (1. - m_omega) * m_effStress;
}

void ElasticDamage::CompEquivanceStrain(double &eps, 
                                        VectorXd &dEpsdStrain,
                                        const std::shared_ptr<Kinematics>&kin)
{
  m_equivanceStrain->CompEquivanceStrain(eps, dEpsdStrain, kin);
}

void ElasticDamage::GetDamage(const double &kappa)
{
  if(kappa <= m_kappa)
  {
    // No Damage
    m_omega = 0.;
    m_dOmegadKappa = 0.;
  }
  else if(kappa < m_kappac)
  {
    // Damage Growth
    double fac = m_kappac / kappa;
    m_omega = fac * (kappa - m_kappa) / (m_kappac - m_kappa);
    m_dOmegadKappa = fac/(m_kappac - m_kappa) - m_omega/kappa;
  }
  else
  {
    // Complete Damage
    m_omega = 1.;
    m_dOmegadKappa = 0.;
  }
}

void ElasticDamage::ComputeDMatrix()
{
  m_D = (1. - m_omega) * m_linearMat->GetTangMatrix();

  if(m_progDam)
    m_D += -m_dOmegadKappa * Math::VecCross(m_effStress, m_dEpsdStrain);
}
