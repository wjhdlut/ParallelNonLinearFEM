#ifndef NONLINEARSOLVER_H
#define NONLINEARSOLVER_H

#include <util/BaseModule.h>

class NonlinearSolver : public BaseModule
{
public:
  NonlinearSolver(const nlohmann::json &props);

  ~NonlinearSolver();

  virtual void Run() override;

private:
  int m_iterMax = 10;
  int m_maxCycle = 20;
  double m_maxLam = 1.0e20;
  double m_tol = 1.0e-3;
};

#endif // NONLINEARSOLVER_H