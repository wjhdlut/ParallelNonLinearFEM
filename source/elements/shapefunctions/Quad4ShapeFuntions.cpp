
#include <elements/shapefunctions/Quad4ShapeFunctions.h>
#include <iostream>

Quad4ShapeFunctions::Quad4ShapeFunctions()
{
  H.resize(4);
  pHpxi.resize(4, 2);
  numOfStress = 3;
}

Quad4ShapeFunctions::~Quad4ShapeFunctions()
{
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