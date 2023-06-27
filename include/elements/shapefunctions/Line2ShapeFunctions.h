#ifndef LINE2SHAPEFUNCTIONS_H
#define LINE2SHAPEFUNCTIONS_H

#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

class Line2ShapeFunctions : public ElementShapeFunctions
{
public:
  Line2ShapeFunctions();
  ~Line2ShapeFunctions();

  virtual void GetShapeFunction(const VectorXd &xi) override;
};
ReflectRegister(Line2ShapeFunctions)

#endif // LINE2SHAPEFUNCTIONS_H