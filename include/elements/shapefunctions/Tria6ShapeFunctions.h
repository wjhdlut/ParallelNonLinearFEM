#ifndef TRIA6SHAPEFUNCTIONS_H
#define TRIA6SHAPEFUNCTIONS_H

#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

/* ----------------------------------------------------
 *  element shape and node order
 *                      3
 *                     / \
 *                    /   \
 *                   /     \
 *                  /       \
 *                 /         \ 
 *                /           \
 *               /             \
 *              1---------------2
 * ----------------------------------------------------- */ 

class Tria6ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct a new Tria 6 Shape Functions object
   * 
   */
  Tria6ShapeFunctions();

  /**
   * @Brief: Destroy the Tria 6 Shape Functions object
   * 
   */
  ~Tria6ShapeFunctions();
  
  /**
   * @Brief: Compute Shape Functions for 8 Nodes Hexahedron
   * 
   * @param xi     [in]  Gauss Points Coordinates
   */
  virtual void GetShapeFunction(const VectorXd &xi) override;
};

ReflectRegister(Tria6ShapeFunctions)

#endif // TRIA6SHAPEFUNCTIONS_H