/**
 * @File Name:     VonMises.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-11-13
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */
#include<iostream>
#include "../../../util/include/Tools.h"
#include "../include/VonMises.h"

VonMises::VonMises(const nlohmann::json &matProps)
{
  Tools::GetParameter(m_nu, "nu", matProps);

  SetAnalyseType(matProps);
}

VonMises::~VonMises()
{
}

void VonMises::SetAnalyseType(const nlohmann::json &matProps)
{
  if(matProps.contains("analyseType"))
  {
    if("PlaneStrain" == matProps.at("analyseType"))
    {
      m_planeStrain = true;
      m_oneVec = VectorXd::Zero(3);
      m_oneVec(0) = 1., m_oneVec(1) = 1.;
    }
    if("PlaneStress" == matProps.at("analyseType"))
    {
      // TO DO IN THE FUTHER
      m_planeStress = true;
    }
    if("AxiSymmetry" == matProps.at("analyseType"))
    {
      m_axiSymmetry = true;
      m_oneVec = VectorXd::Zero(4);
      m_oneVec(0) = 1., m_oneVec(1) = 1., m_oneVec(3) = 1.;
    }
  }
  else
  {
    m_oneVec = VectorXd::Zero(6);
    m_oneVec(0) = 1., m_oneVec(1) = 1., m_oneVec(2) = 1.;
  }
}

double VonMises::CompYieldFunction(const VectorXd &stress)
{
  // std::cout << "stress =\n" << stress << std::endl;
  // std::cout << "m_oneVec =\n " << m_oneVec << std::endl;
  
  if(m_planeStrain)
    return ForPlaneStrain(stress);
  else if(m_planeStress)
    return ForPlaneStress(stress);
  else if(m_axiSymmetry)
    return ForAxiSymmetry(stress);
  else
    return For3D(stress);
  
  return 0.;
}

double VonMises::ForPlaneStrain(const VectorXd &stress)
{
  double stress3 = m_nu * (stress[0] + stress[1]);
  double P = (m_oneVec.dot(stress) + stress3) / 3.;
  VectorXd devStress = stress - P * m_oneVec;
  double devStress3 = stress3 - P; 
  double J2 = sqrt(3./2.*(pow(devStress(0), 2) + pow(devStress(1), 2)
            + pow(devStress3, 2.) + 2. * pow(devStress(2), 2)));
  return J2;
}

double VonMises::ForPlaneStress(const VectorXd &stress)
{
  return 0.;
}

double VonMises::ForAxiSymmetry(const VectorXd &stress)
{
  // std::cout << "stress = " << stress.transpose() << std::endl;
  double P = m_oneVec.dot(stress) / 3.;
  VectorXd devStress = stress - P * m_oneVec;
  // std::cout << "devStress" << devStress.transpose() << std::endl;
  double J2 = sqrt(3./2.*(pow(devStress(0), 2.) + pow(devStress(1), 2.)
              + pow(devStress(3), 2.) + 2.*pow(devStress(2), 2.)));
  return J2;
}

double VonMises::For3D(const VectorXd &stress)
{
  double P = m_oneVec.dot(stress) / 3.;
  VectorXd devStress = stress - P * m_oneVec;
  double J2 =  sqrt(3./2. * (pow(devStress(0), 2) + pow(devStress(1), 2) + pow(devStress(2), 2)
               + 2. * (pow(devStress(3), 2) + pow(devStress(4), 2) + pow(devStress(5), 2))));
  return J2;
}