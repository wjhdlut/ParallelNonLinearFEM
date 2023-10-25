/**
 * @File Name:     Tria3ShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef TRIA3SHAPEFUNCTIONS_H
#define TRIA3SHAPEFUNCTIONS_H

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

class Tria3ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct a new Tria 3 Shape Functions object
   * 
   */
  Tria3ShapeFunctions();

  /**
   * @Brief: Destroy the Tria 3 Shape Functions object
   * 
   */
  ~Tria3ShapeFunctions();
  
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

ReflectRegister(Tria3ShapeFunctions)

#endif // TRIA3SHAPEFUNCTIONS_H