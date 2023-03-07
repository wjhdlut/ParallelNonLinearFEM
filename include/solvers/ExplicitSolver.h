#ifndef EXPLICITSOLVER_H
#define EXPLICITSOLVER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>
#include <petscmat.h>

class ExplicitSolver : public BaseModule
{
public:
  ExplicitSolver(const nlohmann::json &props);

  ~ExplicitSolver();

  virtual void Run() override;

private:
  void PrintStepInfo(const Vec&velo);

private:
  int m_maxCycle;
  double m_dTime;

  std::string m_lam;

  Mat m_lumpedMass;
};
ReflectRegister(ExplicitSolver, const nlohmann::json &)

#endif // EXPLICITSOLVER_H