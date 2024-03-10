/**
 * @File Name:     Math.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <iostream>
#include "../include/Math.h"

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
  if(rowNum*colNum != vec.size()){
    std::cout << "Catch Exception: "
              << "cann't convert vector to matrix"
              << std::endl;
    exit(-1);
  }

  MatrixXd result = MatrixXd::Zero(rowNum, colNum);

  for(int i = 0; i < rowNum; i++)
    for(int j = 0; j < colNum; j++)
      result(i, j) = vec(i*colNum + j);

  return result;
}
}