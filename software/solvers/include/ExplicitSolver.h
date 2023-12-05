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

#include <petscmat.h>

#include "../../util/include/BaseModule.h"
#include "../../util/include/ObjectFactory.h"

/********************************************************************
*   central difference scheme referenced the Box 5.2 in Chapter 5.2  *
*   of the book with name                                            *
*   "Non-Linear Finite Element Analysis of Solids and Structures"    *
*    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel  *
*    John Wiley and Sons, 2012, ISBN 978-0470666449                  *
**********************************************************************/

class ExplicitSolver : public BaseModule
{
public:
  /**
   * @Brief: Construct a new Explicit Solver object
   * 
   * @param props 
   */
  ExplicitSolver(const nlohmann::json &props);
  
  /**
   * @Brief: Destroy the Explicit Solver object
   * 
   */
  ~ExplicitSolver();

  /**
   * @Brief:  Run the Explicit Solver
   * 
   */
  virtual void Run() override;
  
private:
  /**
   * @Brief: Print the Solver Step Information
   * 
   * @param velo 
   */
  void PrintStepInfo(const Vec&velo);

  inline double LamExpression(const double t){
    return 1.e8*(t<1.0e-7);
  }

  void DetermineTimeStepSize();

  void InitialStepComp();

  /**
   * @Brief:         
   * 
   */
  virtual void Initialize(const nlohmann::json &props) override;
  
private:
  int m_maxCycle;
  double m_dTime = 0.;                             // time increment dt^(n-1/2)
  double m_dtScale = 0.9;                          // time scale parameter
  double m_dTime1 = 1.0e6;                         // time increment dt^(n+1/2)
  double m_dTime101d = 1.0e6;
  double m_endTime = 1.0;
  double m_time;
  double m_elemDistortion;

  std::string m_lam;

  Mat m_mass;
  Vec m_lumped;
};
ReflectRegister(ExplicitSolver, const nlohmann::json &)

#endif // EXPLICITSOLVER_H