
#include<util/Math.h>

namespace Math
{
std::vector<std::vector<double>> MatrixAMultB(const std::vector<std::vector<double>> &A,
                                             const std::vector<std::vector<double>> &B)
{
  if(A[0].size() != B.size()) throw "the dimension of A and B is incorrect";
  
  std::vector<std::vector<double>> result;
  std::vector<double> tempVec;
  for(int rowA = 0; rowA < A.size(); rowA++){
    for(int lineB = 0; lineB < B[0].size(); lineB++){
      double temp = 0;
      tempVec.clear();
      for (int rowB = 0; rowB < B.size(); rowB++){
        for (int lineA = 0; lineA < A[0].size(); lineA++){
          temp += A[rowA][lineA] * B[rowB][lineB];
        }
      }
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
    for (int lineB = 0; lineB < B[0].size(); lineB++){
      double temp = 0;
      tempVec.clear();
      for(int rowA = 0; rowA < A.size(); rowA++){
        for(int rowB = 0; rowB< B.size(); rowB++){
          temp += A[rowA][lineA] * B[rowB][lineB];
        }
      }
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
  std::vector<std::vector<double>> invA;
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
}