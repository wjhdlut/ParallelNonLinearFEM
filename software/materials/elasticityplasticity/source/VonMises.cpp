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
  }
  else
  {
    m_oneVec = VectorXd::Zero(6);
    m_oneVec(0) = 1., m_oneVec(1) = 1., m_oneVec(2) = 1.;
  }
}

VonMises::~VonMises()
{
}

double VonMises::CompYieldFunction(const Eigen::VectorXd &stress)
{
  // std::cout << "stress =\n" << stress << std::endl;
  // std::cout << "m_oneVec =\n " << m_oneVec << std::endl;
  
  double stress3 = 0.;
  if(m_planeStrain)
    stress3 = m_nu * (stress[0] + stress[1]);
  
  double P = (m_oneVec.dot(stress) + stress3) / 3.;
  VectorXd devStress = stress - P * m_oneVec;
  
  double devStress3 = 0.;
  if(m_planeStrain)
    devStress3 = stress3 - P;  
  
  double J2 = 0.;
  if(3 == m_oneVec.size())
    J2 =  sqrt(3./2.*(pow(devStress(0), 2) + pow(devStress(1), 2)
        + pow(devStress3, 2.) + 2. * pow(devStress(2), 2)));
  else if (6 == m_oneVec.size())
    J2 =  sqrt(3./2. * (pow(devStress(0), 2) + pow(devStress(1), 2) + pow(devStress(2), 2)
         + 2. * (pow(devStress(3), 2) + pow(devStress(4), 2) + pow(devStress(5), 2))));
  return J2;
}