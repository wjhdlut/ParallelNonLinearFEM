/**
 * @File Name:     RiksSolver.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-22
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef RIKSSOLVER_H
#define RIKSSOLVER_H

#include <petscvec.h>

#include "../../util/include/BaseModule.h"
#include "../../util/include/ObjectFactory.h"

class RiksSolver : public BaseModule
{
public:
  /**
   * @Brief: Construct a new Riks Solver object
   * 
   * @param props 
   */
  RiksSolver(const nlohmann::json &props);

  /**
   * @Brief: Destroy the Riks Solver object
   * 
   */
  ~RiksSolver();

  /**
   * @Brief:  Run the Riks Alc-Length Solver
   * 
   */
  virtual void Run() override;

private:
  void Initialize(const nlohmann::json &props) override;

private:
  int m_optiter = 5;
  int m_iterMax = 10;
  int m_dofCount = 0;
  bool m_fixedStep = false;
  double m_tol = 1.0e-4;
  double m_factor = 1.;
  double m_totalFactor = 1.0;
  double m_maxFactor = 1.0e20;
  double m_maxLam = 1.0e20;
  double m_dLamPrev = 1.;
};

ReflectRegister(RiksSolver, const nlohmann::json &)

#endif // RIKSSOLVER_H