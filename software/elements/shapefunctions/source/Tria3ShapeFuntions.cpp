/**
 * @File Name:     Tria3ShapeFuntions.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/Tria3ShapeFunctions.h"

Tria3ShapeFunctions::Tria3ShapeFunctions()
{
  Initialize();
}

Tria3ShapeFunctions::~Tria3ShapeFunctions()
{}

void Tria3ShapeFunctions::Initialize()
{
  H     = VectorXd::Zero(3);
  pHpxi = MatrixXd::Zero(3, 2);

  numOfStress = 3;
  m_dofType   = {"u", "v"};
}

void Tria3ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  if(2 != xi.size()){
    std::cout << "Catch Exception: "
              << "The isoparamatric coordinate should be 2D for Tria3 element."
              << std::endl;
    exit(-1);
  }
  
  // Calculate shape functions
  H(0) = 1.0-xi(0)-xi(1);
  H(1) = xi(0);
  H(2) = xi(1);

  // Calculate derivatives of shape functions
  pHpxi(0, 0) = -1.0;
  pHpxi(1, 0) =  1.0;
  pHpxi(2, 0) =  0.0;

  pHpxi(0, 1) = -1.0;
  pHpxi(1, 1) =  0.0;
  pHpxi(2, 1) =  1.0;
}