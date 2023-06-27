#include <elements/shapefunctions/Line3ShapeFunctions.h>

Line3ShapeFunctions::Line3ShapeFunctions()
{
  H.resize(3, 0.);
  pHpxi.resize(3, 1);
}

Line3ShapeFunctions::~Line3ShapeFunctions()
{}

void Line3ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  H(0) = 0.5 * (1. - xi(0)) - 0.5 * (1. - xi(0)*xi(0));
  H(1) = 1. - xi(0) * xi(0);
  H(2) = 0.5 * (1. + xi(0)) - 0.5 * (1. - xi(0)*xi(0));

  pHpxi(0, 0) = -0.5 + xi(0);
  pHpxi(1, 0) = -2.  * xi(0);
  pHpxi(2, 0) =  0.5 + xi(0);
}