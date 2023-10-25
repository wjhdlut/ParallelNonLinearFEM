/**
 * @File Name:     LinearSolver.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-22
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>

class LinearSolver : public BaseModule
{
public:
  /**
   * @Brief: Construct a new Linear Solver object
   * 
   * @param props 
   */
  LinearSolver(const nlohmann::json &props);
  
  /**
   * @Brief: Destroy the Linear Solver object
   * 
   */
  ~LinearSolver();

  /**
   * @Brief:  Run the Linear Solver
   * 
   */
  virtual void Run() override;

private:
  int m_iterMax = 10;
  double m_tol = 1.0e-3;
};

ReflectRegister(LinearSolver, const nlohmann::json &)

#endif // LINEARSOLVER_H