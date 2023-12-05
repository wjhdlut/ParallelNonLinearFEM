/**
 * @File Name:     Quad8ShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef QUAD8SHAPEFUNCTION_H
#define QUAD8SHAPEFUNCTION_H

#include "../include/ElementShapeFunctions.h"
#include "../../../util/include/ObjectFactory.h"

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

  /**
   * @Brief: Compute the Boundary Shape Function
   * 
   * @param boundaryH 
   * @param pboundaryHpxi
   */
  virtual void GetBoundaryShapeFunction(VectorXd &boundaryH,
                                        MatrixXd &pboundaryHpxi,
                                        const VectorXd &boundaryXi) override;

  /**
   * @Brief: Get the Boundary Shape Function object
   * 
   * @param boundaryXi 
   * @param boundaryWeight 
   */
  virtual void GetBoundaryIntegrationPoint(MatrixXd &boundaryXi, VectorXd &boundaryWeight) override;

private:
  /**
   * @Brief: Inilize Element Variables
   * 
   */
  virtual void Initialize();

  /**
   * @Brief: Set the Element Node Ordered object
   * 
   * @return std::unordered_map<int, std::vector<int>> 
   */
  virtual std::unordered_map<int, std::vector<int>> SetElemNodeOrdered() override;
};

ReflectRegister(Quad8ShapeFunctions)

#endif // QUAD8SHAPEFUNCTION_H