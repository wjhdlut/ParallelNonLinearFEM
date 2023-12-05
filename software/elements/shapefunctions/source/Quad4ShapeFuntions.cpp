/**
 * @File Name:     Quad4ShapeFuntions.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <iostream>

#include "../include/Quad4ShapeFunctions.h"

Quad4ShapeFunctions::Quad4ShapeFunctions()
{
  Initialize();
}

Quad4ShapeFunctions::~Quad4ShapeFunctions()
{
}

void Quad4ShapeFunctions::Initialize()
{
  H           = VectorXd::Zero(4);
  pHpxi       = MatrixXd::Zero(4, 2);
  numOfStress = 3;
  m_dofType   = {"u", "v"};
}

void Quad4ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  if(2 != xi.size()) throw "The isoparamatric coordinate should be 2D for Quad4 element.";

  // compute shape funtion values at gauss point
  H(0) = 0.25*(1.0-xi(0))*(1.0-xi(1));
  H(1) = 0.25*(1.0+xi(0))*(1.0-xi(1));
  H(2) = 0.25*(1.0+xi(0))*(1.0+xi(1));
  H(3) = 0.25*(1.0-xi(0))*(1.0+xi(1));

  // Calculate derivatives of shape functions 
  pHpxi(0, 0) = -0.25*(1.0-xi(1));
  pHpxi(1, 0) =  0.25*(1.0-xi(1));
  pHpxi(2, 0) =  0.25*(1.0+xi(1));
  pHpxi(3, 0) = -0.25*(1.0+xi(1));

  pHpxi(0, 1) = -0.25*(1.0-xi(0));
  pHpxi(1, 1) = -0.25*(1.0+xi(0));
  pHpxi(2, 1) =  0.25*(1.0+xi(0));
  pHpxi(3, 1) =  0.25*(1.0-xi(0));
}