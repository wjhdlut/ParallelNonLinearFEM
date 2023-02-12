#include <elements/Hexa8ShapeFunctions.h>

Hexa8ShapeFunctions::Hexa8ShapeFunctions()
{
  H.reserve(8);
  
  std::vector<double> temp(3, 0.);
  pHpxi.resize(8, temp);

  numOfStress = 6;
};

Hexa8ShapeFunctions::~Hexa8ShapeFunctions()
{}

void Hexa8ShapeFunctions::GetShapeFunction(const std::vector<double> &xi)
{
  if(3 != xi.size()) throw "The isoparamatric coordinate should be 3D for Hexa8 element.";
  H[0] = 0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0-xi[2]);
  H[1] = 0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[2]);
  H[2] = 0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[2]);
  H[3] = 0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[2]);
  H[4] = 0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[2]);
  H[5] = 0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0+xi[2]);
  H[6] = 0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0+xi[2]);
  H[7] = 0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[2]);

  pHpxi[0][0] =  0.125*(1.0-xi[1])*(1.0-xi[2]);
  pHpxi[1][0] = -0.125*(1.0-xi[1])*(1.0-xi[2]);
  pHpxi[2][0] =  0.125*(1.0+xi[1])*(1.0-xi[2]);
  pHpxi[3][0] = -0.125*(1.0+xi[1])*(1.0-xi[2]);
  pHpxi[4][0] = -0.125*(1.0-xi[1])*(1.0+xi[2]);
  pHpxi[5][0] =  0.125*(1.0-xi[1])*(1.0+xi[2]);
  pHpxi[6][0] =  0.125*(1.0+xi[1])*(1.0+xi[2]);
  pHpxi[7][0] = -0.125*(1.0+xi[1])*(1.0+xi[2]);

  pHpxi[0][1] = -0.125*(1.0-xi[0])*(1.0-xi[2]);
  pHpxi[1][1] = -0.125*(1.0+xi[0])*(1.0-xi[2]);
  pHpxi[2][1] =  0.125*(1.0+xi[0])*(1.0-xi[2]);
  pHpxi[3][1] =  0.125*(1.0-xi[0])*(1.0-xi[2]);
  pHpxi[4][1] = -0.125*(1.0-xi[0])*(1.0+xi[2]);
  pHpxi[5][1] = -0.125*(1.0+xi[0])*(1.0+xi[2]);
  pHpxi[6][1] =  0.125*(1.0+xi[0])*(1.0+xi[2]);
  pHpxi[7][1] =  0.125*(1.0-xi[0])*(1.0+xi[2]);

  pHpxi[0][2] = -0.125*(1.0-xi[0])*(1.0-xi[1]);
  pHpxi[1][2] = -0.125*(1.0+xi[0])*(1.0-xi[1]);
  pHpxi[2][2] = -0.125*(1.0+xi[0])*(1.0+xi[1]);
  pHpxi[3][2] = -0.125*(1.0-xi[0])*(1.0+xi[1]);
  pHpxi[4][2] =  0.125*(1.0-xi[0])*(1.0-xi[1]);
  pHpxi[5][2] =  0.125*(1.0+xi[0])*(1.0-xi[1]);
  pHpxi[6][2] =  0.125*(1.0+xi[0])*(1.0+xi[1]);
  pHpxi[7][2] =  0.125*(1.0-xi[0])*(1.0+xi[1]);
}