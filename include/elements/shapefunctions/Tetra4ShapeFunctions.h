/**
 * @File Name:     Tetra4ShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef TETRA4SHAPEFUNCTIONS_H
#define TETRA4SHAPEFUNCTIONS_H

#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

class Tetra4ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct a new Tetra 4 Shape Functions object
   * 
   */
  Tetra4ShapeFunctions();

  /**
   * @Brief: Destroy the Tetra 4 Shape Functions object
   * 
   */
  ~Tetra4ShapeFunctions();
  
  /**
   * @Brief:  Compute the Shape Functions
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

ReflectRegister(Tetra4ShapeFunctions)

#endif // TETRA4SHAPEFUNCTIONS_H