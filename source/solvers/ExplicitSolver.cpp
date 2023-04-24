#include <solvers/ExplicitSolver.h>
#include <util/DataStructure.h>
#include <iostream>
#include <iomanip>

ExplicitSolver::ExplicitSolver(const nlohmann::json &props) : BaseModule(props)
{
  Tools::GetParameter(m_maxCycle, "maxCycle", m_myProps);
  Tools::GetParameter(m_dTime, "dtime", m_myProps);
  Tools::GetParameter(m_lam, "lam", m_myProps);
  
  VecDuplicate(GlobalData::GetInstance()->m_state, &m_lumped);
  GlobalData::GetInstance()->m_elements->AssembleMassMatrix(m_mass, m_lumped);

  Tools::PrintVecIntoFile(m_lumped, "lumped_mass.txt");
}

ExplicitSolver::~ExplicitSolver()
{
  VecDestroy(&m_lumped);
  MatDestroy(&m_mass);
}

void ExplicitSolver::Run()
{
  DetermineTimeStepSize();
  GlobalData::GetInstance()->m_cycle += 1;
  GlobalData::GetInstance()->m_time += m_dTime;

  // Get the Displacement, Velocity, Acceleration
  Vec &disp = GlobalData::GetInstance()->m_state;
  Vec &velo = GlobalData::GetInstance()->m_velo;
  Vec &acce = GlobalData::GetInstance()->m_acce;
  Vec &dDisp = GlobalData::GetInstance()->m_Dstate;
  
  // Get the Internal Force and External Force
  Vec &fint = GlobalData::GetInstance()->m_fint;
  Vec &fhat = GlobalData::GetInstance()->m_fhat;
 
 // update velocity, acceleration, displacement and increment displacement
  VecAXPY(velo, 0.25*(m_dTime+m_dTime1), acce);
  VecAXPY(disp, m_dTime, velo);
  VecCopy(velo, dDisp);
  VecScale(dDisp, m_dTime);
  

  GlobalData::GetInstance()->m_elements->AssembleInternalForce(fint);

  double lam = LamExpression(GlobalData::GetInstance()->m_time);
  GlobalData::GetInstance()->m_dofs->SetConstrainFactor(lam);

  Vec res;
  VecDuplicate(fhat, &res);
  VecCopy(fhat, res);
  VecScale(res, lam);
  VecAXPY(res, -1., fint);
  // std::cout << "res = " << std::endl;
  // VecView(res, PETSC_VIEWER_STDOUT_WORLD);

  GlobalData::GetInstance()->m_dofs->Solve(m_lumped, res, acce);
  // std::cout << "acce = " << std::endl;
  // VecView(acce, PETSC_VIEWER_STDOUT_WORLD);

  VecAXPY(velo, 0.25*(m_dTime + m_dTime1), acce);
  // std::cout << "velo = " << std::endl;
  // VecView(velo, PETSC_VIEWER_STDOUT_WORLD);

  m_dTime = m_dTime1;

  GlobalData::GetInstance()->m_elements->CommitHistory();

  PrintStepInfo(velo);

  if(GlobalData::GetInstance()->m_cycle == m_maxCycle)
    GlobalData::GetInstance()->m_active = false;

  VecDestroy(&res);
}

void ExplicitSolver::PrintStepInfo(const Vec&velo)
{
  if(0 == GlobalData::GetInstance()->m_cycle % 20 
  || 1 == GlobalData::GetInstance()->m_cycle)
  {
    std::cout << "  Cycle     Time         Kin.Energy" << std::endl;
    std::cout << "  ---------------------------------------" << std::endl;
  }

  double energy;
  Vec tempVec;
  VecDuplicate(velo, &tempVec);
  VecPointwiseMult(tempVec, m_lumped, velo);
  VecDot(tempVec, velo, &energy);
  
  std::cout << std::setw(7) << GlobalData::GetInstance()->m_cycle;

  std::cout.precision(3);
  std::cout << std::setw(13) << std::setiosflags(std::ios::scientific)
            << GlobalData::GetInstance()->m_time
            << std::setw(13) << std::setiosflags(std::ios::scientific) << 0.5 * energy << std::endl;
}

void ExplicitSolver::DetermineTimeStepSize(){
    GlobalData::GetInstance()->m_elements->ReturnCurentTimeIncrementPara(m_dTime, 
                                                     m_dTime1, m_elemDistortion);
    if (abs(m_dTime - 1.0e-6) > 1.0e-5){
      m_dTime = std::min(m_dTime * m_dtScale, 1.05 * m_dTime1);
      m_dTime1 = m_dTime;
    }
  }