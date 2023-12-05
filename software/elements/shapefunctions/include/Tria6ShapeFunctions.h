/**
 * @File Name:     Tria6ShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef TRIA6SHAPEFUNCTIONS_H
#define TRIA6SHAPEFUNCTIONS_H

#include "../include/ElementShapeFunctions.h"
#include "../../../util/include/ObjectFactory.h"

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

private:
  /**
   * @Brief: Inilize Element Variables
   * 
   */
  virtual void Initialize();
};

ReflectRegister(Tria6ShapeFunctions)

#endif // TRIA6SHAPEFUNCTIONS_H