/**
 * @File Name:     Line2ShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef LINE2SHAPEFUNCTIONS_H
#define LINE2SHAPEFUNCTIONS_H

#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

/* ----------------------------------------------------
 *  element shape and node order
 *              1----------------------2
 * ----------------------------------------------------- */ 

class Line2ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct a new Line 2 Shape Functions object
   * 
   */
  Line2ShapeFunctions();

  /**
   * @Brief: Destroy the Line 2 Shape Functions object
   * 
   */
  ~Line2ShapeFunctions();
  
  /**
   * @Brief: Compute the Shape Functions
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
ReflectRegister(Line2ShapeFunctions)

#endif // LINE2SHAPEFUNCTIONS_H