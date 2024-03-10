/**
 * @File Name:     VonMisesEquivantStrain.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2024-01-23
 * 
 * @Copyright Copyright (c) 2024 JianHuaWang
 * 
 */
#include <iostream>

#include "../include/VonMisesEquivantStrain.h"
#include "../../../util/include/Tools.h"

VonMisesEquivantStrain::VonMisesEquivantStrain(const nlohmann::json &prop)
                       : EquivantStrainBase(prop)
{
  // Read Some Basic Variables
  ReadParameters();

  // Initialize ...
  Initialize();
}

VonMisesEquivantStrain::~VonMisesEquivantStrain()
{}

void VonMisesEquivantStrain::CompEquivanceStrain(double &eps,
                                                 VectorXd &dEpsdStrain,
                                                 const std::shared_ptr<Kinematics>&kin)
{
  if(m_planeStrainFlag)
    dEpsdStrain = PlaneStrainEquivanceStrain(eps, kin);
  else if(m_planeStressFlag)
    dEpsdStrain = PlaneStressEquivandeStrain(eps, kin);
  else
    dEpsdStrain = ThreeDimensionEquivandeStrain(eps, kin);
}

void VonMisesEquivantStrain::ReadParameters()
{
  Tools::GetParameter(m_k, "k", m_props);
  Tools::GetParameter(m_nu, "nu", m_props);
}

Vector3d VonMisesEquivantStrain::PlaneStrainEquivanceStrain(double &eps,
                                                        const std::shared_ptr<Kinematics>&kin)
{
  double ezz = m_nu / (m_nu - 1.) * (kin->strain(0) + kin->strain(1));

  // the First Invariant of Strain Tensor
  double I1 = kin->strain(0) + kin->strain(1) + ezz;

  // the second Invariant of Deviatoric Strain Tensor
  double temp1 = pow(kin->strain(0), 2.) + pow(kin->strain(1), 2.) + pow(ezz, 2.);
  double temp2 = kin->strain(0) * kin->strain(1) + kin->strain(0) * ezz + kin->strain(1) * ezz;
  double J2 =  ( temp1 - temp2) / 3. + pow(kin->strain(2), 2.);

  // Compute the derivative of equivalent strain w.r.t. strain array
  Vector3d dExxdStrain = Vector3d(1, 0, 0);
  Vector3d dEyydStrain = Vector3d(0, 1, 0);
  Vector3d dExydStrain = 0.5 * Vector3d(0, 0, 1);
  Vector3d dEzzdStrain = m_c * (dExxdStrain + dEyydStrain);

  // Compute the derivative of strain first invariant w.r.t. strain array
  Vector3d dI1dStrain  = dExxdStrain + dEyydStrain + dEzzdStrain;
  
  // Compute the derivative of deviatoric strain second invariant w.r.t. strain array
  Vector3d dJ2dStrain = Vector3d::Zero();
  dJ2dStrain = m_sc * (2.*kin->strain(0)-kin->strain(1)-ezz) * dExxdStrain;
  dJ2dStrain += m_sc * (2.*kin->strain(1)-kin->strain(0)-ezz) * dEyydStrain;
  dJ2dStrain += m_sc * (2.*ezz-kin->strain(0)-kin->strain(1)) * dEzzdStrain;
  dJ2dStrain += 2. * kin->strain(2) * dExydStrain;

  double disc = pow(m_a2 * I1, 2.) + m_a3*J2;
  double dDiscdI1 = 2. * pow(m_a2, 2.) * I1;
  double dDiscdJ2 = m_a3;

  if(disc < 1.e-16)
  {
    double temp3 = 0.;
    Vector3d dTemp3dStrain = Vector3d(m_a4, m_a4, 2.*m_a3);

    // Equivalent Strain
    eps = m_a1 * (m_a2*I1 + temp3);
    
    return m_a1 * m_a2 * dI1dStrain + m_a1 * dTemp3dStrain;
  }
  else
  {
    double temp3 = sqrt(disc);
    double dTemp3dI1 = 0.5/temp3*dDiscdI1;
    double dTemp3dJ2 = 0.5/temp3*dDiscdJ2;

    double dEpsdI1 = m_a1 * (m_a2 + dTemp3dI1);
    double dEpsdJ2 = m_a1 * dTemp3dJ2;
    
    // Equivalent Strain
    eps = m_a1 * (m_a2 * I1 + temp3);

    return dEpsdI1 * dI1dStrain + dEpsdJ2 * dJ2dStrain;
  }
}

Vector3d VonMisesEquivantStrain::PlaneStressEquivandeStrain(double &eps,
                                                            const std::shared_ptr<Kinematics>&kin)
{
  // To Do Future
  return Vector3d::Zero();
}

VectorXd VonMisesEquivantStrain::ThreeDimensionEquivandeStrain(double &eps,
                                                               const std::shared_ptr<Kinematics>&kin)
{
// To Do Future
return VectorXd::Zero(6);
}

void VonMisesEquivantStrain::Initialize()
{
  m_c  = m_nu / (m_nu - 1.);
  m_a1 = 1. / (2. * m_k);
  m_a2 = (m_k - 1.) / (1. - 2. * m_nu);
  m_a3 = 12. * m_k / pow(1. + m_nu, 2.);
  m_a4 = sqrt(pow(m_a2, 2.) + m_a3*m_sc);
}