#ifndef QUAD8SHAPEFUNCTION_H
#define QUAD8SHAPEFUNCTION_H

#include <elements//shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

/* ----------------------------------------------------
 *  element shape and node order
 *              7----------6-----------5
 *              |                      |
 *              |                      |
 *              |                      |
 *              8                      4
 *              |                      |
 *              |                      |
 *              |                      |
 *              1----------2-----------3
 * ----------------------------------------------------- */ 

class Quad8ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct Shape Functions object for 8 Nodes Quadrangle
   * 
   */
  Quad8ShapeFunctions();

  /**
   * @Brief: Destroy Shape Functions object for 8 Nodes Quadrangle
   * 
   */
  ~Quad8ShapeFunctions();

  /**
   * @Brief:  Compute the Shape Functions
   * 
   * @param xi   [in]  Gauss Points Coordinates
   */
  virtual void GetShapeFunction(const VectorXd &xi) override;
};

ReflectRegister(Quad8ShapeFunctions)

#endif // QUAD8SHAPEFUNCTION_H