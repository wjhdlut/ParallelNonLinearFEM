
#include <util/Math.h>
#include <iostream>

namespace Math
{
MatrixXd VecCross(const VectorXd &vecA, const VectorXd &vecB)
{
  MatrixXd result = MatrixXd::Zero(vecA.size(), vecB.size());
  for(int i = 0; i < vecA.size(); i++)
    for(int j = 0; j < vecB.size(); j++)
      result(i, j) = vecA(i) * vecB(j);

  return result;
}

MatrixXd ConvertVecToMat(const int rowNum, const int colNum, const VectorXd &vec)
{
  if(rowNum*colNum != vec.size()) throw "cann't convert vector to matrix";

  MatrixXd result = MatrixXd::Zero(rowNum, colNum);

  for(int i = 0; i < rowNum; i++)
    for(int j = 0; j < colNum; j++)
      result(i, j) = vec(i*colNum + j);

  return result;
}
}