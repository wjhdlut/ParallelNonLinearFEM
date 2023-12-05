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

#include "../include/VonMises.h"

VonMises::VonMises(const nlohmann::json &matProps)
{
  if(matProps.contains("analyseType"))
  {
    if("PlaneStrain" == matProps.at("analyseType"))
    {
      m_oneVec = VectorXd::Zero(3);
      m_oneVec(0) = 1., m_oneVec(1) = 1.;
    }
    if("PlaneStress" == matProps.at("analyseType"))
    {
      // TO DO IN THE FUTHER
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
  double P = m_oneVec.dot(stress) / 3.;
  VectorXd devStress = stress - P * m_oneVec;
  
  double J2 = 0.;
  if(3 == m_oneVec.size())
    J2 =  sqrt(3./2.*(pow(devStress(0), 2) + pow(devStress(1), 2)
        + 2. * pow(devStress(2), 2)));
  else if (6 == m_oneVec.size())
    J2 =  sqrt(3./2. * (pow(devStress(0), 2) + pow(devStress(1), 2) + pow(devStress(2), 2)
         + 2. * (pow(devStress(3), 2) + pow(devStress(4), 2) + pow(devStress(5), 2))));
  return J2;
}