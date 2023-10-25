/**
 * @File Name:     Line3ShapeFunctions.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <elements/shapefunctions/Line3ShapeFunctions.h>

Line3ShapeFunctions::Line3ShapeFunctions()
{
  Initialize();
}

Line3ShapeFunctions::~Line3ShapeFunctions()
{}

void Line3ShapeFunctions::Initialize()
{
  H         = VectorXd::Zero(3);
  pHpxi     = MatrixXd::Zero(3, 1);
  m_dofType = {"u", "v"};
}

void Line3ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  H(0) = 0.5 * (1. - xi(0)) - 0.5 * (1. - xi(0)*xi(0));
  H(1) = 1. - xi(0) * xi(0);
  H(2) = 0.5 * (1. + xi(0)) - 0.5 * (1. - xi(0)*xi(0));

  pHpxi(0, 0) = -0.5 + xi(0);
  pHpxi(1, 0) = -2.  * xi(0);
  pHpxi(2, 0) =  0.5 + xi(0);
}