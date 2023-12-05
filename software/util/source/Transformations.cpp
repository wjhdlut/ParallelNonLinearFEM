/**
 * @File Name:     Transformations.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/Transformations.h"
#include "../include/Math.h"

namespace Transformations
{
VectorXd ToElementCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords)
{
  return VecToElementCoordinates(a, elemNodeCoords);
}

MatrixXd ToElementCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords)
{
  return MatToElementCoordinates(A, elemNodeCoords);
}

VectorXd ToGlobalCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords)
{
  return VecToGLobalCoordinates(a, elemNodeCoords);
}

MatrixXd ToGlobalCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords)
{
  return MatToGlobalCoordinates(A, elemNodeCoords);
}

VectorXd VecToElementCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords)
{
  MatrixXd R = GetRotationMatrix(elemNodeCoords);

  VectorXd aBar = VectorXd::Zero(a.size());

  if(0 != a.size()%R.rows()) throw "Vector does not have the right shape to be rotated";

  VectorXd tempVec = VectorXd::Zero(R.rows());
  for(int i = 0; i < a.size()/R.rows(); i++)
  {
    for(int j = 0; j < R.rows(); j++)
      tempVec(j) = a(i*R.rows() + j);
    
    tempVec = R * tempVec;
    
    for(int j = 0; j < R.rows(); j++)
      aBar(i*R.rows() + j) = tempVec(j);
  }

  return aBar;
}

MatrixXd MatToElementCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords)
{
  MatrixXd R = GetRotationMatrix(elemNodeCoords);

  VectorXd tempVec = VectorXd::Zero(A.cols());
  MatrixXd aBar(A.rows(), A.cols());

  MatrixXd tempMat = MatrixXd::Zero(R.rows(), R.rows());
  
  if(0 != A.cols() % R.rows() || 0 != A.rows() % R.rows())
    throw "Matrix does not have the right shape to be rotated";
  
  for(int i = 0; i < A.rows()/R.rows(); i++)
  {
    for(int j = 0; j < A.cols()/R.rows(); j++)
    {
      for(int ii = 0; ii < R.rows(); ii++)
        for(int jj = 0; jj < R.rows(); jj++)
          tempMat(ii, jj) = A(i*R.rows() + ii, j*R.rows() + jj);

      tempMat = R * tempMat * R.transpose();

      for(int ii = 0; ii < R.rows(); ii++)
        for(int jj = 0; jj < R.rows(); jj++)
          aBar(i*R.rows() + ii, j*R.rows() + jj) = tempMat(ii, jj);
    }
  }

  return aBar;
}

MatrixXd GetRotationMatrix(const MatrixXd &elemNodeCoords)
{
  if(0 == elemNodeCoords.rows()) return MatrixXd::Zero(0, 0);

  if(2 != elemNodeCoords.cols()) throw "Rotation matrix only implemented for 2D situation";


  // Compute the undeformed element length
  double elemLength = (elemNodeCoords.row(1) - elemNodeCoords.row(0)).norm();

  // rotate a globdal coordinate to an element coordinate
  double sin_alpha = (elemNodeCoords(1, 1) - elemNodeCoords(0, 1))/elemLength;
  double cos_alpha = (elemNodeCoords(1, 0) - elemNodeCoords(0, 0))/elemLength;
  
  Matrix2d A;
  A << cos_alpha, sin_alpha, -sin_alpha, cos_alpha;
  return A;
}

VectorXd VecToGLobalCoordinates(const VectorXd &a, const MatrixXd &elemNodeCoords)
{
  MatrixXd R = GetRotationMatrix(elemNodeCoords);

  VectorXd aBar = VectorXd::Zero(a.size());

  if(0 != a.size()%R.rows()) throw "Vector does not have the right shape to be rotated";

  VectorXd tempVec = VectorXd::Zero(R.rows());

  for(int i = 0; i < a.size()/R.rows(); i++)
  {
    for(int j = 0; j < R.rows(); j++)
      tempVec(j) = a(i*R.rows() + j);

    tempVec = R * tempVec;

    for(int j = 0; j < R.rows(); j++)
      aBar(i*R.rows() + j) = tempVec(j);
  }

  return aBar;
}

MatrixXd MatToGlobalCoordinates(const MatrixXd &A, const MatrixXd &elemNodeCoords)
{
  MatrixXd R = GetRotationMatrix(elemNodeCoords);

  MatrixXd aBar = MatrixXd::Zero(A.rows(), A.cols());
  MatrixXd tempMat = MatrixXd::Zero(R.rows(), R.cols());
  
  if(0 != A.cols() % R.rows() || 0 != A.rows() % R.cols())
    throw "Matrix does not have the right shape to be rotated";
  
  for(int i = 0; i < A.cols()/R.rows(); i++)
  {
    for(int j = 0; j < A.rows()/R.rows(); j++)
    {

      for(int ii = 0; ii < R.rows(); ii++)
        for(int jj = 0; jj < R.rows(); jj++)
          tempMat(ii, jj) = A(i*R.rows() + ii, j*R.rows() + jj);

      tempMat = R.transpose() * tempMat * R;

      for(int ii = 0; ii < R.rows(); ii++)
        for(int jj = 0; jj < R.rows(); jj++)
          aBar(i*R.rows() + ii, j*R.rows() + jj) = tempMat(ii, jj);
    }
  }

  return aBar;
}
}