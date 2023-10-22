#include <elements/shapefunctions/Tetra4ShapeFunctions.h>

Tetra4ShapeFunctions::Tetra4ShapeFunctions()
{
  Initialize();
}

Tetra4ShapeFunctions::~Tetra4ShapeFunctions()
{}

void Tetra4ShapeFunctions::Initialize()
{
  H     = VectorXd::Zero(4);
  pHpxi = MatrixXd::Zero(4, 3);

  numOfStress = 6;
  m_dofType   = {"u", "v"};
}

void Tetra4ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  if(3 != xi.size()) throw "The isoparamatric coordinate should be 3D for Hexa8 element.";

  // Calculate shape functions
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

  pHpxi(0, 2) = -0.25*(1.0-xi(0));
  pHpxi(1, 2) = -0.25*(1.0+xi(0));
  pHpxi(2, 2) =  0.25*(1.0+xi(0));
  pHpxi(3, 2) =  0.25*(1.0-xi(0));
}
