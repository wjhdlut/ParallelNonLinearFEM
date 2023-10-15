#include <elements/shapefunctions/Quad8ShapeFunctions.h>

Quad8ShapeFunctions::Quad8ShapeFunctions()
{
  Initialize();
};

Quad8ShapeFunctions::~Quad8ShapeFunctions()
{}

void Quad8ShapeFunctions::Initialize()
{
  H.resize(8);
  pHpxi.resize(8, 2);
  numOfStress = 3;
  // dofType = {"u", "v"};
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