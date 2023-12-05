#ifndef NONLINEARSOLVER_H
#define NONLINEARSOLVER_H

#include "../../util/include/BaseModule.h"
#include "../../util/include/ObjectFactory.h"

class NonlinearSolver : public BaseModule
{
public:
  /**
   * @Brief: Construct a new Nonlinear Solver object
   * 
   * @param props 
   */
  NonlinearSolver(const nlohmann::json &props);

  /**
   * @Brief: Destroy the Nonlinear Solver object
   * 
   */
  ~NonlinearSolver();

  /**
   * @Brief: Nonlinear Solver Algothrim
   * 
   */
  virtual void Run() override;

  /**
   * @Brief:  Read Solver Data form Input file
   * 
   * @param fileName 
   */
  virtual void ReadData(const std::string &fileName) override;

private:
  /**
   * @Brief: Initialize Some Basic Variables
   * 
   * @param props 
   */
  virtual void Initialize(const nlohmann::json &props) override;

  /**
   * @Brief: Update the Load Factor Lambda
   * 
   * @param numOfCycle 
   */
  void UpdateLoadFactor(const int &numOfCycle);

private:
  // int m_iterMax = 10;
  int m_maxCycle = 0;
  double m_maxLam = 1.0e20;
  // double m_tol = 1.0e-3;

  std::vector<int> m_iterMaxVec;
  std::vector<double> m_lamVec;
  std::vector<double> m_tolVec;
};

ReflectRegister(NonlinearSolver, const nlohmann::json&)

#endif // NONLINEARSOLVER_H