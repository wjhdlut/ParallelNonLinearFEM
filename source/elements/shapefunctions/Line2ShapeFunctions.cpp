#include <elements/shapefunctions/Line2ShapeFunctions.h>

Line2ShapeFunctions::Line2ShapeFunctions()
{
  H.resize(2);

  std::vector<double> temp(1, 0.);
  pHpxi.resize(2, temp);
}

Line2ShapeFunctions::~Line2ShapeFunctions()
{}

void Line2ShapeFunctions::GetShapeFunction(const std::vector<double> &xi)
{
  H[0] = 0.5 * (1. - xi[0]);
  H[1] = 0.5 * (1. + xi[0]);

  pHpxi[0][0] = -0.5;
  pHpxi[1][0] = 0.5;
}