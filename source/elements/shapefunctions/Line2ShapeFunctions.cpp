#include <elements/shapefunctions/Line2ShapeFunctions.h>

Line2ShapeFunctions::Line2ShapeFunctions()
{
  H.resize(2);
  pHpxi.resize(2, 1);
}

Line2ShapeFunctions::~Line2ShapeFunctions()
{}

void Line2ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  H(0) = 0.5 * (1. - xi(0));
  H(1) = 0.5 * (1. + xi(0));

  pHpxi(0, 0) = -0.5;
  pHpxi(1, 0) = 0.5;
}