#include <iostream>
#include <solvers/NonlinearSolver.h>
#include <util/DataStructure.h>
#include <petscksp.h>

NonlinearSolver::NonlinearSolver(const nlohmann::json &props)
                 : BaseModule(props)
{
  GlobalData::GetInstance()->m_lam = 0.;
}

NonlinearSolver::~NonlinearSolver()
{}

void NonlinearSolver::Run()
{
  GlobalData::GetInstance()->m_cycle += 1;
  GlobalData::GetInstance()->m_lam = 1.0 * GlobalData::GetInstance()->m_cycle;

  Vec &a = GlobalData::GetInstance()->m_state;
  Vec &Da = GlobalData::GetInstance()->m_Dstate;
  Vec &fhat = GlobalData::GetInstance()->m_fhat;

  Vec &fint = GlobalData::GetInstance()->m_fint;

  Vec fext;
  VecDuplicate(fint, &fext);

  std::cout << "=================================" << std::endl;
  std::cout << " Load step " << GlobalData::GetInstance()->m_cycle << std::endl;
  std::cout << "=================================" << std::endl;
  std::cout << "  NR iter : L2-norm residual" << std::endl;

  GlobalData::GetInstance()->m_iiter = 0;

  Mat K;
  GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);

  double error = 1.;

  while(error > m_tol)
  {
    GlobalData::GetInstance()->m_iiter += 1;

    GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);
  }
}