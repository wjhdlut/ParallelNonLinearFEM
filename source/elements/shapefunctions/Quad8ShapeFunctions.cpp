/**
 * @File Name:     Quad8ShapeFunctions.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <elements/shapefunctions/Quad8ShapeFunctions.h>
#include <elements/shapefunctions/Line3ShapeFunctions.h>
#include <util/ShapeFunctions.h>

Quad8ShapeFunctions::Quad8ShapeFunctions()
{
  Initialize();
};

Quad8ShapeFunctions::~Quad8ShapeFunctions()
{}

void Quad8ShapeFunctions::Initialize()
{
  H = VectorXd::Zero(8);
  pHpxi = MatrixXd::Zero(8, 2);
  numOfStress = 3;
  m_dofType = {"u", "v"};
}

void Quad8ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  if(2 != xi.size()) throw "The isoparamatric coordinate should be 2D for Quad8 element.";

  // compute shape funtion values at gauss point
  H(0) = -0.25*(1.0-xi(0))*(1.0-xi(1))*(1.0+xi(0)+xi(1));
  H(1) =  0.5 *(1.0-xi(0))*(1.0+xi(0))*(1.0-xi(1));
  H(2) = -0.25*(1.0+xi(0))*(1.0-xi(1))*(1.0-xi(0)+xi(1));
  H(3) =  0.5 *(1.0+xi(0))*(1.0+xi(1))*(1.0-xi(1));
  H(4) = -0.25*(1.0+xi(0))*(1.0+xi(1))*(1.0-xi(0)-xi(1));
  H(5) =  0.5 *(1.0-xi(0))*(1.0+xi(0))*(1.0+xi(1));
  H(6) = -0.25*(1.0-xi(0))*(1.0+xi(1))*(1.0+xi(0)-xi(1));
  H(7) =  0.5 *(1.0-xi(0))*(1.0+xi(1))*(1.0-xi(1));

  // Calculate derivatives of shape functions 
  pHpxi(0, 0) = -0.25*(-1.0+xi(1))*( 2.0*xi(0)+xi(1));
  pHpxi(1, 0) =  xi(0)*(-1.0+xi(1));
  pHpxi(2, 0) =  0.25*(-1.0+xi(1))*(-2.0*xi(0)+xi(1));
  pHpxi(3, 0) = -0.5 *(1.0+xi(1))*(-1.0+xi(1));
  pHpxi(4, 0) =  0.25*( 1.0+xi(1))*( 2.0*xi(0)+xi(1));
  pHpxi(5, 0) = -xi(0)*(1.0+xi(1));
  pHpxi(6, 0) = -0.25*( 1.0+xi(1))*(-2.0*xi(0)+xi(1));
  pHpxi(7, 0) = 0.5*(1.0+xi(1))*(-1.0+xi(1));

  pHpxi(0, 1) = -0.25*(-1.0+xi(0))*( xi(0)+2.0*xi(1));
  pHpxi(1, 1) =  0.5 *( 1.0+xi(0))*(-1.0+xi(0));
  pHpxi(2, 1) =  0.25*( 1.0+xi(0))*(-xi(0)+2.0*xi(1));
  pHpxi(3, 1) = -xi(1)*(1.0+xi(0));
  pHpxi(4, 1) =  0.25*( 1.0+xi(0))*( xi(0)+2.0*xi(1));
  pHpxi(5, 1) = -0.5 *( 1.0+xi(0))*(-1.0+xi(0));
  pHpxi(6, 1) = -0.25*(-1.0+xi(0))*(-xi(0)+2.0*xi(1));
  pHpxi(7, 1) =  xi(1)*(-1.0+xi(0));
}

void Quad8ShapeFunctions::GetBoundaryShapeFunction(VectorXd &boundaryH, 
                                                   MatrixXd &pboundaryHpxi,
                                                   const VectorXd &boundaryXi)
{
  Line3ShapeFunctions *res = new Line3ShapeFunctions();
  res->GetShapeFunction(boundaryXi);
  boundaryH = res->H;
  pboundaryHpxi = res->pHpxi;
  delete res;
}

std::unordered_map<int, std::vector<int>> Quad8ShapeFunctions::SetElemNodeOrdered()
{
  std::unordered_map<int, std::vector<int>> elemNodeOrdered;
  elemNodeOrdered.insert(std::pair<int, std::vector<int>>(1, {1, 2, 3, 0, 0, 0, 0, 0}));
  elemNodeOrdered.insert(std::pair<int, std::vector<int>>(2, {0, 0, 1, 2, 3, 0, 0, 0}));
  elemNodeOrdered.insert(std::pair<int, std::vector<int>>(3, {0, 0, 0, 0, 1, 2, 3, 0}));
  elemNodeOrdered.insert(std::pair<int, std::vector<int>>(4, {3, 0, 0, 0, 0, 0, 1, 2}));

  return elemNodeOrdered;
}

void Quad8ShapeFunctions::GetBoundaryIntegrationPoint(MatrixXd &boundaryXi, VectorXd &boundaryWeight)
{
  ShapeFunctions::GetIntegrationPoints(boundaryXi, boundaryWeight, "Line3", -1);
}