#ifndef HEXA8SHAPEFUNCTION_H
#define HEXA8SHAPEFUNCTION_H

#include <elements/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

class Hexa8ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct Shape Functions object for 8 Nodes Hexahedron
   * 
   */
  Hexa8ShapeFunctions();

  /**
   * @Brief: Destroy Shape Functions object for 8 Nodes Hexahedron
   * 
   */
  ~Hexa8ShapeFunctions();

  /**
   * @Brief:  Compute Shape Functions for 8 Nodes Hexahedron
   * 
   * @param xi     [in]  Gauss Points Coordinates
   */
  virtual void GetShapeFunction(const std::vector<double> &xi) override;
};

ReflectRegister(Hexa8ShapeFunctions);

#endif // HEXA8SHAPEFUNCTION_H