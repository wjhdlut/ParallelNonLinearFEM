
#include <util/Math.h>
#include <iostream>
#include <math.h>
#include <numeric>

namespace Math
{
std::vector<std::vector<double>> MatrixAMultB(const std::vector<std::vector<double>> &A,
                                             const std::vector<std::vector<double>> &B)
{
  if(A[0].size() != B.size()) throw "the dimension of A and B is incorrect";
  
  std::vector<std::vector<double>> result;
  std::vector<double> tempVec;
  for(int rowA = 0; rowA < A.size(); rowA++){
    tempVec.clear();
    for(int lineB = 0; lineB < B[0].size(); lineB++){
      double temp = 0;
      for (int rowB = 0; rowB < B.size(); rowB++)
        temp += A[rowA][rowB] * B[rowB][lineB];
      tempVec.emplace_back(temp);
    }
    result.emplace_back(tempVec);
  }

  return result;
}

std::vector<std::vector<double>> MatrixATransMultB(const std::vector<std::vector<double>> &A,
                                                   const std::vector<std::vector<double>> &B)
{
  if(A.size() != B.size()) throw "the dimension of A and B is incorrect";

  std::vector<std::vector<double>> result;
  std::vector<double> tempVec;
  for(int lineA = 0; lineA < A[0].size(); lineA++){
    tempVec.clear();
    for (int lineB = 0; lineB < B[0].size(); lineB++){
      double temp = 0;
      for(int row = 0; row < A.size(); row++)
        temp += A[row][lineA] * B[row][lineB];
      tempVec.emplace_back(temp);
    }
    result.emplace_back(tempVec);
  }
  return result;
}

std::vector<std::vector<double>> MatrixAMultBTrans(const std::vector<std::vector<double>> &A,
                                                   const std::vector<std::vector<double>> &B)
{
  if(A[0].size() != B[0].size()) throw "the dimension of A and B is incorrect";

  std::vector<std::vector<double>> result;
  std::vector<double> tempVec;

  for(int rowA = 0; rowA < A.size(); rowA++)
  {
    tempVec.clear();
    for(int rowB = 0; rowB < B.size(); rowB++)
    {
      double temp = 0;
      for(int line = 0; line < A[0].size(); line++)
        temp += A[rowA][line] * B[rowB][line];
      tempVec.emplace_back(temp);
    }
    result.emplace_back(tempVec);
  }

  return result;
}

double MatrixDet(const std::vector<std::vector<double>> &A)
{
  if(A.size() != A[0].size()) throw "Array to compute det is not square matrix";
  if(A.size() > 3) throw "inability to compute det of matrix with rank more than 3";

  double det = 0;
  if(A.size() == 2)
    return A[1][1]*A[0][0] - A[0][1]*A[1][0];

  if(A.size() == 3){
    double temp = A[0][0]*A[1][1]*A[2][2] + A[1][0]*A[2][1]*A[0][2] + A[2][0]*A[0][1]*A[1][2]
                - A[2][0]*A[1][1]*A[0][2] - A[0][1]*A[1][0]*A[2][2] - A[0][0]*A[2][1]*A[1][2];
    return temp;
  }
  return 0;
}

std::vector<std::vector<double>> MatrixInverse(const std::vector<std::vector<double>>&A)
{
  if(A.size() != A[0].size()) throw "Array for inverse is not square matrix";
  if(A.size() > 3) throw "inability to compute inverse of matrix with rank more than 3";

  double detA = MatrixDet(A);
  std::vector<std::vector<double>> invA = A;
  if(2 == A.size()){
    invA[0][0] = A[1][1]/detA;
    invA[1][1] = A[0][0]/detA;

    invA[0][1] = -A[0][1]/detA;
    invA[1][0] = -A[1][0]/detA;
  }

  if(3 == A.size()){
    invA[0][0] = ( A[1][1]*A[2][2] - A[1][2]*A[2][1])/detA;
    invA[0][1] = (-A[1][0]*A[2][2] + A[1][2]*A[2][0])/detA;
    invA[0][2] = ( A[1][0]*A[2][1] - A[1][1]*A[2][0])/detA;

    invA[1][0] = (-A[0][1]*A[2][2] + A[0][2]*A[2][1])/detA;
    invA[1][1] = ( A[0][0]*A[2][2] + A[0][2]*A[2][0])/detA;
    invA[1][2] = (-A[0][0]*A[2][1] + A[0][1]*A[2][0])/detA;

    invA[2][0] = ( A[0][1]*A[1][2] - A[1][1]*A[0][2])/detA;
    invA[2][1] = (-A[0][0]*A[1][2] + A[1][0]*A[0][2])/detA;
    invA[2][2] = ( A[0][0]*A[1][1] - A[0][1]*A[1][0])/detA;
  }
  
  return invA;
}

std::vector<std::vector<double>> MatrixEye(const int dim)
{
  std::vector<double> rowData(dim, 0);
  std::vector<std::vector<double>> eye(dim, rowData);
  for(int i = 0; i < dim; i++)
    eye[i][i] = 1.;
  return eye;
}

std::vector<std::vector<double>> MatrixAdd(const double scale,
                                           const std::vector<std::vector<double>>&A,
                                           const std::vector<std::vector<double>>&B)
{
  if(A.size() != B.size() || A[0].size() != B[0].size())
    throw "the size of matrix A and B should be same for add";
  
  std::vector<double> temp(A[0].size(), 0);
  std::vector<std::vector<double>> result(A.size(), temp);
  for (int i = 0; i < A.size(); i++)
    for(int j = 0; j < A[0].size(); j++)
      result[i][j] = A[i][j] + scale * B[i][j];
  return result;
}

std::vector<double> MatrixAMultVecB(const std::vector<std::vector<double>>&A,
                                     const std::vector<double>&b)
{
  if(A[0].size() != b.size()) throw "Matrix and Vector is incompatible";
  std::vector<double> c(A.size(), 0.);

  for(int i = 0; i < A.size(); i++)
    for(int j = 0; j < b.size(); j++)
      c[i] += A[i][j] * b[j];
  
  return c;
}

Matrix MatrixZeros(const int rowNum, const int lineNum)
{
  std::vector<double> temp(lineNum, 0.);
  Matrix A(rowNum, temp);

  return A;
}

Matrix MatrixScale(const double scale, const Matrix&A)
{
  Matrix B = A;
  for(int i = 0; i < A.size(); i++)
    for(int j = 0; j < A[0].size(); j++)
      B[i][j] = scale * A[i][j];
  
  return B;
}

std::vector<double> MatrixATransMultVecB(const Matrix&A,
                                          const std::vector<double>&b)
{
  if(A.size() != b.size()) throw "Matrix and Vector is incompatible";
  std::vector<double> c(A[0].size(), 0.);
  for(int i = 0; i < A[0].size(); i++)
    for(int j = 0; j < A.size(); j++)
      c[i] += A[j][i] * b[j];

  return c;
}

std::vector<double> VecScale(const double scale, const std::vector<double> &A)
{
  std::vector<double> b = A;
  for(int i = 0; i < A.size(); i++)
    b[i] = scale * A[i];

  return b;
}

std::vector<double> VecAdd(const double scale, const std::vector<double>&a, const std::vector<double>&b)
{
  if(a.size() != b.size()) throw "size of two vectors to add should be same.";
  std::vector<double> c = a;

  for(int i = 0; i < a.size(); i++)
    c[i] = a[i] + scale * b[i];
  
  return c;
}

std::vector<double> ConvertMatrixToVec(const std::vector<std::vector<double>>&A)
{
  std::vector<double> B;
  for(auto iRowValue : A)
    B.insert(B.end(), iRowValue.begin(), iRowValue.end());
  
  return B;
}

Matrix VecOuter(const std::vector<double>&A, const std::vector<double>&B)
{
  std::vector<double> temp(B.size(), 0.);
  Matrix C(A.size(), temp);
  for(int row = 0; row < A.size(); row++)
    for(int line = 0; line < B.size(); line++)
      C[row][line] = A[row] * B[line];

  return C;
}

void MatrixOutput(const std::vector<std::vector<double>> &A)
{
  for(int i = 0; i < A.size(); i++)
  {
    for(int j = 0; j < A[0].size(); j++)
    {
      std::cout << A[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}

double VecNorm(const std::vector<double> &A)
{
  double temp = 0;
  for(auto a : A)
    temp += a*a;
  return sqrt(temp);
}

double VecSum(const std::vector<double>&a)
{
  return std::accumulate(a.begin(), a.end(), 0.);
} 
}