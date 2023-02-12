#ifndef NONLINEARSOLVER_H
#define NONLINEARSOLVER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>

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

private:
  int m_iterMax = 10;
  int m_maxCycle = 20;
  double m_maxLam = 1.0e20;
  double m_tol = 1.0e-3;
};

ReflectRegister(NonlinearSolver, const nlohmann::json&)

#endif // NONLINEARSOLVER_H