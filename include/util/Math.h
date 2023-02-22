#ifndef MATH_H
#define MATH_H

#include <vector>

typedef std::vector<std::vector<double>> Matrix;

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

/**
 * @Brief:  form unit matrix with rank dim
 * 
 * @param dim 
 * @return std::vector<std::vector<double>> 
 */
std::vector<std::vector<double>> MatrixEye(const int dim);

/**
 * @Brief: compute C = A + scale * B;
 * 
 * @param scale 
 * @param A 
 * @param B 
 * @return std::vector<std::vector<double>> 
 */
std::vector<std::vector<double>> MatrixAdd(const double scale,
                                 const std::vector<std::vector<double>>&A,
                                 const std::vector<std::vector<double>>&B);

std::vector<double> VecAdd(const double scale, const std::vector<double>&a, const std::vector<double>&b);

/**
 * @Brief:   compute result A*b
 * 
 * @param A   Matrix
 * @param b   vector
 * @return std::vector<double> 
 */
std::vector<double> MatrixAMultVecB(const std::vector<std::vector<double>>&A,
                                    const std::vector<double>&b);

/**
 * @Brief:   compute result A^T*b
 * 
 * @param A   Matrix
 * @param b   vector
 * @return std::vector<double> 
 */
std::vector<double> MatrixTAMultVecB(const Matrix&A,
                                     const std::vector<double>&b);
/**
 * @Brief:  form a zero matrix 
 * 
 * @param rowNum 
 * @param lineNum 
 * @return Matrix 
 */
Matrix MatrixZeros(const int rowNum, const int lineNum);

/**
 * @Brief:  compute result B = scale * A for matrix
 * 
 * @param scale 
 * @param A 
 * @return Matrix 
 */
Matrix MatrixScale(const double scale, const Matrix&A);

/**
 * @Brief:  compute result B = scale * A for vector
 * 
 * @param scale 
 * @param A 
 * @return std::vector<double> 
 */
std::vector<double> VecScale(const double scale, const std::vector<double> &A);

/**
 * @Brief: Compute the outer product of two vectors.
 * 
 * @param A = [A0, A1, ..., Am]
 * @param B = [B0, B1, ..., Bn]
 * @return Matrix = [A0*B0, A0*B1, A0*B2, ..., A0*Bn;
 *                   A1*B0, A1*B1, A1*B2, ..., A1*Bn;
 *                   ...;
 *                   Am*B0, Am*B1, Am*B2, ..., Am*Bn;]
 */
Matrix VecOuter(const std::vector<double>&A, const std::vector<double>&B);

/**
 * @Brief:  Convert Matrix to Vector
 * 
 * @param A 
 * @return std::vector<double> 
 */
std::vector<double> ConvertMatrixToVec(const std::vector<std::vector<double>> &A);

/**
 * @Brief: Output Matrix
 * 
 * @param A 
 */
void MatrixOutput(const std::vector<std::vector<double>>&A);
};

#endif // MATH_H
