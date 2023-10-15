
#include <elements/shapefunctions/Tria3ShapeFunctions.h>

Tria3ShapeFunctions::Tria3ShapeFunctions()
{
  Initialize();
}

Tria3ShapeFunctions::~Tria3ShapeFunctions()
{}

void Tria3ShapeFunctions::Initialize()
{
  H.resize(3);
  pHpxi.resize(3, 2);

  numOfStress = 3;
}

void Tria3ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  if(2 != xi.size()) throw "The isoparamatric coordinate should be 2D for Tria3 element.";
  
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