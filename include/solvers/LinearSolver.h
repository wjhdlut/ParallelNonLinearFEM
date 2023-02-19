#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>

class LinearSolver : public BaseModule
{
public:
  LinearSolver(const nlohmann::json &props);

  ~LinearSolver();

  virtual void Run() override;

private:
  int m_iterMax = 10;
  double m_tol = 1.0e-3;
};

ReflectRegister(LinearSolver, const nlohmann::json &)

#endif // LINEARSOLVER_H