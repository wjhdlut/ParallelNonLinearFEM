/**
 * @File Name:     ExplicitSolver.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-03-13
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef EXPLICITSOLVER_H
#define EXPLICITSOLVER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>
#include <petscmat.h>

/* central difference scheme referenced the Box 5.2 in Chapter 5.2
*  of the book with name
*   "Non-Linear Finite Element Analysis of Solids and Structures"    *
*    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel  *
*    John Wiley and Sons, 2012, ISBN 978-0470666449                  *
*/

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