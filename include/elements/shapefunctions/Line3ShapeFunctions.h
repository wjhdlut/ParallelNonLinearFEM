#ifndef LINE3SHAPEFUNCTIONS_H
#define LINE3SHAPEFUNCTIONS_H

#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

class Line3ShapeFunctions : public ElementShapeFunctions
{
public:
  Line3ShapeFunctions();
  ~Line3ShapeFunctions();

  virtual void GetShapeFunction(const VectorXd &xi);
};
ReflectRegister(Line3ShapeFunctions)

#endif // LINE3SHAPEFUNCTIONS_H