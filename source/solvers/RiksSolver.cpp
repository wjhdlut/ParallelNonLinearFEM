#include <solvers/RiksSolver.h>
#include <util/DataStructure.h>
#include <iostream>

RiksSolver::RiksSolver(const nlohmann::json &props) : BaseModule(props)
{
  GlobalData::GetInstance()->m_lam = 1.0;
}

RiksSolver::~RiksSolver()
{

}

void RiksSolver::Run()
{
  GlobalData::GetInstance()->m_cycle += 1;

  Vec &a = GlobalData::GetInstance()->m_state;
  Vec &Da = GlobalData::GetInstance()->m_Dstate;
  Vec &fhat = GlobalData::GetInstance()->m_fhat;
  Vec &fint = GlobalData::GetInstance()->m_fint;

  Vec fext, da1, daPrev, res;
  VecDuplicate(fhat, &fext);
  VecCopy(fhat, fext);
  VecScale(fext, GlobalData::GetInstance()->m_lam);

  VecDuplicate(fhat, &da1);
  VecDuplicate(fhat, &daPrev);
  VecDuplicate(fhat, &res);

  std::cout << "======================================" << std::endl;
  std::cout << " Load step " << GlobalData::GetInstance()->m_cycle << std::endl;
  std::cout << "======================================" << std::endl;
  std::cout << "  iter # : L2-norm residual" << std::endl;

  // Initialize Newton-Raphson Iteration Parameters
  double error = 1.;
  GlobalData::GetInstance()->m_iiter = 0;

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  // Predictor
  Mat K;
  double dLam = 0;
  double dLamPrev = 1.0;
  if(1 == GlobalData::GetInstance()->m_cycle)
  {
    GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);
    GlobalData::GetInstance()->m_dofs->Solve(K, fext, da1, ksp);
    dLam = GlobalData::GetInstance()->m_lam;
  }
  else
  {
    VecCopy(da1, daPrev);
    VecScale(da1, m_factor);

    dLam = m_factor * dLamPrev;
    GlobalData::GetInstance()->m_lam += dLam;
  }

  VecAXPY(a, 1., da1);
  VecAXPY(Da, 1., da1);
  
  GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);

  VecWAXPY(res, -1., fext, fint);

  Vec tempVecD1, tempVecD2, dDa;
  VecDuplicate(fext, &tempVecD1);
  VecDuplicate(fext, &tempVecD2);
  VecDuplicate(fext, &dDa);

  double dDLam = 0., tempDValue1= 0., tempDValue2 = 0.;
  while(error > m_tol)
  {
    GlobalData::GetInstance()->m_iiter += 1;
    
    GlobalData::GetInstance()->m_dofs->Solve(K, fext, tempVecD1, ksp);
    GlobalData::GetInstance()->m_dofs->Solve(K, res, tempVecD2, ksp);

    VecDot(da1, tempVecD2, &tempDValue2);
    VecDot(da1, tempVecD1, &tempDValue1);
    dDLam = - tempDValue2 / tempDValue1;

    VecWAXPY(dDa, dDLam, tempVecD1, tempVecD2);

    dLam += dDLam;
    GlobalData::GetInstance()->m_lam += dDLam;

    VecAXPY(a, 1., dDa);
    VecAXPY(Da, 1., dDa);

    GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);

    VecCopy(fhat, fext);
    VecScale(fext, GlobalData::GetInstance()->m_lam);
    VecWAXPY(res, -1., fext, fint);

    GlobalData::GetInstance()->m_dofs->Norm(res, tempDValue1);
    GlobalData::GetInstance()->m_dofs->Norm(fext, tempDValue2);
    error = tempDValue1 / tempDValue2;

    // Print Iterative Information
    std::cout << "  Iter " << GlobalData::GetInstance()->m_iiter
              << " : error = " << error << std::endl;

    if(m_iterMax == GlobalData::GetInstance()->m_iiter)
      throw "Newton-Raphson iterations did not converge!";
  }
  // Print Converged Information
  std::cout << "--------------------------------------" << std::endl;
  std::cout << " Converged in " << GlobalData::GetInstance()->m_iiter << " iterations" << std::endl;

  GlobalData::GetInstance()->m_elements->CommitHistory();

  VecCopy(fint, GlobalData::GetInstance()->m_fint);

  if(!m_fixedStep)
  {
    m_factor = pow(0.5, 0.25*(GlobalData::GetInstance()->m_iiter - m_optiter));
    m_totalFactor *= m_factor;
  }

  if(m_totalFactor > m_maxFactor)
    m_factor = 1.;

  VecCopy(Da, daPrev);
  dLamPrev = dLam;

  if(GlobalData::GetInstance()->m_lam > m_maxLam ||
     GlobalData::GetInstance()->m_cycle > 1000)
     GlobalData::GetInstance()->m_active = false;
}