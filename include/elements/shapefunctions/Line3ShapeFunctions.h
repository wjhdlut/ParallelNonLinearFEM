/**
 * @File Name:     Line3ShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef LINE3SHAPEFUNCTIONS_H
#define LINE3SHAPEFUNCTIONS_H

#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

/* ----------------------------------------------------
 *  element shape and node order
 *              1-----------2-----------3
 * ----------------------------------------------------- */ 

class Line3ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct a new Line 3 Shape Functions object
   * 
   */
  Line3ShapeFunctions();

  /**
   * @Brief: Destroy the Line 3 Shape Functions object
   * 
   */
  ~Line3ShapeFunctions();
  
  /**
   * @Brief: Compute the Shape Functions
   * 
   * @param xi   [in]  Gauss Points Coordinates
   */
  virtual void GetShapeFunction(const VectorXd &xi);

private:
  /**
   * @Brief: Inilize Element Variables
   * 
   */
  virtual void Initialize();
};
ReflectRegister(Line3ShapeFunctions)

#endif // LINE3SHAPEFUNCTIONS_H