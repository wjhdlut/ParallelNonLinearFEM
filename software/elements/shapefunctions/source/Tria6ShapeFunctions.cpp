/**
 * @File Name:     Tria6ShapeFunctions.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/Tria6ShapeFunctions.h"

Tria6ShapeFunctions::Tria6ShapeFunctions()
{
  Initialize();
}

Tria6ShapeFunctions::~Tria6ShapeFunctions()
{}

void Tria6ShapeFunctions::Initialize()
{
  H     = VectorXd::Zero(6);
  pHpxi = MatrixXd::Zero(6, 2);

  numOfStress = 3;
  m_dofType   = {"u", "v"};
}

void Tria6ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  if(2 != xi.size()) throw "The isoparamatric coordinate should be 2D for Tria6 element.";

  // Calculate shape functions
  H(0) = 1.0-xi[0]-xi[1]-2.0*xi[0]*(1.0-xi[0]-xi[1])-2.0*xi[1]*(1.0-xi[0]-xi[1]);
  H(1) = xi[0]-2.0*xi[0]*(1.0-xi[0]-xi[1])-2.0*xi[0]*xi[1];
  H(2) = xi[1]-2.0*xi[0]*xi[1]-2.0*xi[1]*(1.0-xi[0]-xi[1]);
  H(3) = 4.0*xi[0]*(1.0-xi[0]-xi[1]);
  H(4) = 4.0*xi[0]*xi[1];
  H(5) = 4.0*xi[1]*(1.0-xi[0]-xi[1]);

  // Calculate derivatives of shape functions
  pHpxi(0, 0) = -1.0-2.0*(1.0-xi[0]-xi[1])+2.0*xi[0]+2.0*xi[1];
  pHpxi(1, 0) =  1.0-2.0*(1.0-xi[0]-xi[1])+2.0*xi[0]-2.0*xi[1];
  pHpxi(2, 0) =  0.0;
  pHpxi(3, 0) =  4.0*(1.0-xi[0]-xi[1])-4.0*xi[0];
  pHpxi(4, 0) =  4.0*xi[1];
  pHpxi(5, 0) = -4.0*xi[1];

  pHpxi(0, 1) = -1.0+2.0*xi[0]-2.0*(1.0-xi[0]-xi[1])+2.0*xi[1];
  pHpxi(1, 1) =  0.0;
  pHpxi(2, 1) =  1.0-2.0*xi[0]-2.0*(1.0-xi[0]-xi[1])+2.0*xi[1];
  pHpxi(3, 1) = -4.0*xi[0];
  pHpxi(4, 1) =  4.0*xi[0];
  pHpxi(5, 1) =  4.0*(1.0-xi[0]-xi[1])-4.0*xi[1];
}