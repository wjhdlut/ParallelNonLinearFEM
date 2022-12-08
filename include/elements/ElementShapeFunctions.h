
#ifndef ELEMENTSHAPEFUNCTIONS_H
#define ELEMENTSHAPEFUNCTIONS_H

#include <vector>

class ElementShapeFunctions
{
public:
  ElementShapeFunctions(){}
  virtual ~ElementShapeFunctions(){}

  virtual void GetShapeFunction(const std::vector<double> &xi) = 0;

public:
  std::vector<double> H;
  std::vector<std::vector<double>> pHpxi;
};

#endif // ELEMENTSHAPEFUNCTIONS_H
