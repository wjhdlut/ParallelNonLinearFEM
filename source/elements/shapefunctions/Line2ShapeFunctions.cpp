#include <elements/shapefunctions/Line2ShapeFunctions.h>

Line2ShapeFunctions::Line2ShapeFunctions()
{
  Initialize();
}

Line2ShapeFunctions::~Line2ShapeFunctions()
{}

void Line2ShapeFunctions::Initialize()
{
  H         = VectorXd::Zero(2);
  pHpxi     = MatrixXd::Zero(2, 1);
  m_dofType = {"u", "v"};
}

void Line2ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  H(0) = 0.5 * (1. - xi(0));
  H(1) = 0.5 * (1. + xi(0));

  pHpxi(0, 0) = -0.5;
  pHpxi(1, 0) = 0.5;
}

