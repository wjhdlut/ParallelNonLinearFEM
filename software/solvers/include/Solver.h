/**
 * @File Name:     Solver.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-22
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef SOLVERS_H
#define SOLVERS_H

#include <memory>

#include "../../util/include/BaseModule.h"

class Solver
{
public:
  /**
   * @Brief:  Construct a new Solver object
   * 
   */
  Solver();

  /**
   * @Brief:  Destroy the Solver object
   * 
   */
  virtual ~Solver();
  
  /**
   * @Brief:  Solvering Process
   * 
   */
  void Run();

private:
  void Initialize();

protected:
  std::shared_ptr<BaseModule> m_solver;
};



#endif // SOLVERS_H