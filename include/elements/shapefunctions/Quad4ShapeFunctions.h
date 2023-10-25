/**
 * @File Name:     Quad4ShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef QUAD4SHAPEFUNCTIONS_H
#define QUAD4SHAPEFUNCTIONS_H

#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

/* ----------------------------------------------------
 *  element shape and node order
 *              4----------------------3
 *              |                      |
 *              |                      |
 *              |                      |
 *              |                      |
 *              |                      |
 *              |                      |
 *              |                      |
 *              1----------------------2
 * ----------------------------------------------------- */ 

class Quad4ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct a new Quad 4 Shape Functions object
   * 
   */
  Quad4ShapeFunctions();

  /**
   * @Brief: Destroy the Quad 4 Shape Functions object
   * 
   */
  ~Quad4ShapeFunctions();

  /**
   * @Brief:  Compute the Shape Functions
   * 
   * @param xi   [in]  Gauss Points Coordinates
   */
  virtual void GetShapeFunction(const VectorXd &xi) override;

private:
  /**
   * @Brief: Inilize Element Variables
   * 
   */
  virtual void Initialize();
};

ReflectRegister(Quad4ShapeFunctions)

#endif // QUAD4SHAPEFUNCTIONS_H