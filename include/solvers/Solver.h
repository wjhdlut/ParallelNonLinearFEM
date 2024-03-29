

#ifndef SOLVERS_H
#define SOLVERS_H

#include <memory>
#include <util/BaseModule.h>

class Solver
{
public:
  Solver();
  virtual ~Solver();
  
  void Run();

protected:
  std::shared_ptr<BaseModule> m_solver;
};



#endif // SOLVERS_H