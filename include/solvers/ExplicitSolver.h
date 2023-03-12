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

  inline double LamExpression(const double t){
    return 1.e8*(t<1.0e-7);
  }

private:
  int m_maxCycle;
  double m_dTime;

  std::string m_lam;

  Mat m_mass;
  Vec m_lumped;
};
ReflectRegister(ExplicitSolver, const nlohmann::json &)

#endif // EXPLICITSOLVER_H