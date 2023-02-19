
#include <elements/shapefunctions/Tria3ShapeFunctions.h>

Tria3ShapeFunctions::Tria3ShapeFunctions()
{
  H.reserve(3);
  std::vector<double> temp(2);
  pHpxi.resize(3, temp);

  numOfStress = 3;
}

Tria3ShapeFunctions::~Tria3ShapeFunctions()
{}

void Tria3ShapeFunctions::GetShapeFunction(const std::vector<double> &xi)
{
  if(2 != xi.size()) throw "The isoparamatric coordinate should be 2D for Tria3 element.";
  
  // Calculate shape functions
  H[0] = 1.0-xi[0]-xi[1];
  H[1] = xi[0];
  H[2] = xi[1];

  // Calculate derivatives of shape functions
  pHpxi[0][0] = -1.0;
  pHpxi[1][0] =  1.0;
  pHpxi[2][0] =  0.0;

  pHpxi[0][1] = -1.0;
  pHpxi[1][1] =  0.0;
  pHpxi[2][1] =  1.0;
}