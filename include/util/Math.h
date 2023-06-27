#ifndef MATH_H
#define MATH_H

#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace Math
{
/**
 * @Brief: Compute the cross of two vectors
 * 
 * @param a 
 * @param b 
 * @return MatrixXd 
 */
MatrixXd VecCross(const VectorXd &vecA, const VectorXd &VecB);

/**
 * @Brief: Convert vector to matrix
 * 
 * @param rowNum 
 * @param colNum 
 * @param vec 
 * @return MatrixXd 
 */
MatrixXd ConvertVecToMat(const int rowNum, const int colNum, const VectorXd &vec);
};

#endif // MATH_H
