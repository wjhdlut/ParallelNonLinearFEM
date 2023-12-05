/**
 * @File Name:     Transformations.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace Transformations
{
/**
 * @Brief: Transform vector to element coornidate system
 * 
 * @param a 
 * @param elemNodeCoords 
 * @return VectorXd 
 */
VectorXd ToElementCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords);

/**
 * @Brief: Transform matrix to element coordinate system
 * 
 * @param A 
 * @param elemNodeCoords 
 * @return MatrixXd 
 */
MatrixXd ToElementCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords);

/**
 * @Brief: Transform vector to global coordinate system
 * 
 * @param a 
 * @param elemNodeCoords 
 * @return VectorXd 
 */
VectorXd ToGlobalCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords);

/**
 * @Brief: Transform matrix to global coordinate system
 * 
 * @param A 
 * @param elemNodeCoords 
 * @return MatrixXd 
 */
MatrixXd ToGlobalCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords);

/**
 * @Brief:  Vector Variables Transform to Element Coordinate System
 * 
 * @param a 
 * @param elemNodeCoords 
 * @return VectorXd 
 */
VectorXd VecToElementCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords);

/**
 * @Brief:  Matrix Variables Transform to Element Coordinate System
 * 
 * @param A 
 * @param elemNodeCoords 
 * @return MatrixXd 
 */
MatrixXd MatToElementCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords);

/**
 * @Brief: Vector Variables Transform to Global Coordinate System
 * 
 * @param a 
 * @param elemNodeCoords 
 * @return VectorXd 
 */
VectorXd VecToGLobalCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords);

/**
 * @Brief: Matrix Variables Transform to Global Coordinate System
 * 
 * @param A 
 * @param elemNodeCoords 
 * @return MatrixXd 
 */
MatrixXd MatToGlobalCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords);

/**
 * @Brief:  Compute the Rotation Matrix
 * 
 * @param elemNodeCoords 
 * @return MatrixXd 
 */
MatrixXd GetRotationMatrix(const MatrixXd &elemNodeCoords);

}

#endif // TRANSFORMATIONS_H