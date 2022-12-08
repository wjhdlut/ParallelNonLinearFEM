#ifndef MATH_H
#define MATH_H

#include <vector>

namespace Math
{
/**
 * @Brief: compute matrix mlutiplication A*B
 * 
 * @param A 
 * @param B 
 * @return std::vector<std::vector<double>> 
 */
std::vector<std::vector<double>> MatrixAMultB(const std::vector<std::vector<double>> &A,
                                             const std::vector<std::vector<double>> &B);

/**
 * @Brief: compute matrix mlutiplication A^T*B
 * 
 * @param A 
 * @param B 
 * @return std::vector<std::vector<double>> 
 */
std::vector<std::vector<double>> MatrixATransMultB(const std::vector<std::vector<double>> &A,
                                                  const std::vector<std::vector<double>> &B);
/**
 * @Brief:  compute the det of square matrix
 * 
 * @param A 
 * @return double 
 */
double MatrixDet(const std::vector<std::vector<double>> &A);

/**
 * @Brief:  compute the inverse of square matrix
 * 
 * @param A 
 * @return std::vector<std::vector<double>> 
 */
std::vector<std::vector<double>> MatrixInverse(const std::vector<std::vector<double>> &A);
};

#endif // MATH_H
