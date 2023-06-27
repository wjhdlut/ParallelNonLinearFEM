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
 * @Brief:         
 * 
 * @param a 
 * @param elemNodeCoords 
 * @return VectorXd 
 */
VectorXd VecToElementCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords);

MatrixXd MatToElementCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords);

VectorXd VecToGLobalCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords);

MatrixXd MatToGlobalCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords);

MatrixXd GetRotationMatrix(const MatrixXd &elemNodeCoords);

}

#endif // TRANSFORMATIONS_H