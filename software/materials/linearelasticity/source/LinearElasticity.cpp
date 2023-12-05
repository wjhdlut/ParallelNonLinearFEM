/**
 * @File Name:     LinearElasticity.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/LinearElasticity.h"

LinearElasticity::LinearElasticity(const nlohmann::json &props) : BaseMaterial(props)
{
  m_E = SetMaterialParamter("E");
  m_nu = SetMaterialParamter("nu");
  ComputeDMatrix();
}

LinearElasticity::~LinearElasticity()
{}

void LinearElasticity::ComputeDMatrix()
{
  if (m_props.contains("analyseType"))
  {
    if ("PlaneStrain" == m_props.at("analyseType"))
      ForPlaneStrain();
    if ("PlaneStress" == m_props.at("analyseType"))
      ForPlaneStress();
  }
  else{
    For3D();
  }
}

void LinearElasticity::For3D()
{
  double fac = 1.0 / (2.0 * m_nu * m_nu + m_nu - 1.0 );
  m_D = MatrixXd::Zero(6, 6);

  m_D(0, 0) = fac * m_E * ( m_nu - 1.0 );
  m_D(0, 1) = -1.0 * fac * m_E * m_nu;
  m_D(0, 2) = m_D(0, 1);
  m_D(1, 0) = m_D(0, 1);
  m_D(1, 1) = m_D(0, 0);
  m_D(1, 2) = m_D(0, 1);
  m_D(2, 0) = m_D(0, 1);
  m_D(2, 1) = m_D(0, 1);
  m_D(2, 2) = m_D(0, 0);
  m_D(3, 3) = m_E / ( 2.0 + 2.0 * m_nu );
  m_D(4, 4) = m_D(3, 3);
  m_D(5, 5) = m_D(3, 3);
}

void LinearElasticity::ForPlaneStrain()
{
  m_D = MatrixXd::Zero(3, 3);
  m_D(0, 0) = m_E*(1.-m_nu)/((1+m_nu)*(1.-2.*m_nu));
  m_D(0, 1) = m_D(0, 0)*m_nu/(1-m_nu);
  m_D(1, 0) = m_D(0, 1);
  m_D(1, 1) = m_D(0, 0);
  m_D(2, 2) = m_D(0, 0)*0.5*(1.-2.*m_nu)/(1.-m_nu);
}

void LinearElasticity::ForPlaneStress()
{
  m_D = MatrixXd::Zero(3, 3);
  m_D(0, 0) = m_E/(1. - m_nu * m_nu);
  m_D(0, 1) = m_D(0, 0)*m_nu;
  m_D(1, 0) = m_D(0, 1);
  m_D(1, 1) = m_D(0, 0);
  m_D(2, 2) = m_E/2./(1. + m_nu);
}