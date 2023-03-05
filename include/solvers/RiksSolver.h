#ifndef RIKSSOLVER_H
#define RIKSSOLVER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>
#include <petscvec.h>

class RiksSolver : public BaseModule
{
public:
  RiksSolver(const nlohmann::json &props);

  ~RiksSolver();

  virtual void Run() override;

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