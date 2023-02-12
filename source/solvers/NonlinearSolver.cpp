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
  VecDuplicate(fhat, &fext);
  VecCopy(fhat, fext);
  VecScale(fext, GlobalData::GetInstance()->m_lam);

  std::cout << "=================================" << std::endl;
  std::cout << " Load step " << GlobalData::GetInstance()->m_cycle << std::endl;
  std::cout << "=================================" << std::endl;
  std::cout << "  NR iter : L2-norm residual" << std::endl;

  GlobalData::GetInstance()->m_iiter = 0;

  Mat K;
  GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);
  // MatView(K, PETSC_VIEWER_STDOUT_WORLD);

  Vec dF, da;
  VecDuplicate(fext, &dF);
  VecWAXPY(dF, -1., fint, fext);
  VecDuplicate(fext, &da);

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  double norm = 0., error = 1.;
  while(error > m_tol)
  {
    GlobalData::GetInstance()->m_iiter += 1;

    GlobalData::GetInstance()->m_dofs->Solve(K, dF, da, ksp);
    // std::cout << "da = " << std::endl;
    // VecView(da, PETSC_VIEWER_STDOUT_WORLD);

    // update the increment displacement and total displacement
    VecAXPY(Da, 1., da);
    VecAXPY(a, 1., da);

    GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);
    VecWAXPY(dF, -1., fint, fext);
    // std::cout << "fint = " << std::endl;
    // VecView(fint, PETSC_VIEWER_STDOUT_WORLD);

    VecNorm(fext, NORM_2, &norm);

    if(norm < 1.0e-16){
      GlobalData::GetInstance()->m_dofs->Norm(dF, error);
    }
    else{
      GlobalData::GetInstance()->m_dofs->Norm(dF, error);

      error /= norm;
    }

    std::cout << "  Iter " << GlobalData::GetInstance()->m_iiter
              << " : error = " << error << std::endl;
    
    if(GlobalData::GetInstance()->m_iiter == m_iterMax)
      throw "Newton-Raphson iterations did not converge!";
  }

  GlobalData::GetInstance()->m_elements->CommitHistory();

  VecSet(Da, 0.);

  if(GlobalData::GetInstance()->m_cycle == m_maxCycle 
  || GlobalData::GetInstance()->m_lam > m_maxLam)
    GlobalData::GetInstance()->m_active = false;

  VecDestroy(&da);
  VecDestroy(&dF);
  VecDestroy(&fext);
  MatDestroy(&K);
  KSPDestroy(&ksp);
}