#ifndef HEXA8SHAPEFUNCTION_H
#define HEXA8SHAPEFUNCTION_H

#include <elements/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

class Hexa8ShapeFunctions : public ElementShapeFunctions
{
public:
  Hexa8ShapeFunctions();
  ~Hexa8ShapeFunctions();

  virtual void GetShapeFunction(const std::vector<double> &xi) override;
};

ReflectRegister(Hexa8ShapeFunctions);

#endif // HEXA8SHAPEFUNCTION_H