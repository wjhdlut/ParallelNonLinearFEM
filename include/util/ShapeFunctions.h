/**
 * @File Name:     ShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef SHAPEFUNCTIONS_H
#define SHAPEFUNCTIONS_H

#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace ShapeFunctions
{
/**
 * @Brief: Get Integration Point Coordinate and Weight Bassed
 *         on Gauss Integration strategy 
 * 
 * @param xi 
 * @param weight 
 * @param order 
 */
void GaussScheme(VectorXd&xi, VectorXd&weight, const int order);

/**
 * @Brief: Get Integration Point Coordinate and Weight Bassed
 *         on Tria Integration strategy
 * 
 * @param order 
 */
void TriaScheme(const int order);

/**
 * @Brief: Get the Element Type Based on Nodal Dimension and 
 *         Number of Node in Element
 * 
 * @param elemCoords 
 * @return std::string 
 */
std::string GetElemType(const MatrixXd&elemCoords);

/**
 * @Brief:  Get the Integration Points For Various Element
 * 
 * @param xi 
 * @param weight 
 * @param elemType 
 * @param order 
 */
void GetIntegrationPoints(MatrixXd&xi, VectorXd&weight, const std::string&elemType,
                          const int order = 0);
}

#endif // SHAPEFUNCTIONS_H
