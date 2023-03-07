#include <solvers/ExplicitSolver.h>
#include <util/DataStructure.h>
#include <petscksp.h>
#include <iomanip>

ExplicitSolver::ExplicitSolver(const nlohmann::json &props) : BaseModule(props)
{
  if(props["solver"].contains("maxCycle"))
    m_maxCycle = props["solver"].at("maxCycle");
  
  if(props["solver"].contains("dtime"))
    m_dTime = props["solver"].at("dtime");

  if(props["solver"].contains("lam"))
    m_lam = props["solver"].at("lam");
}

ExplicitSolver::~ExplicitSolver()
{}

void ExplicitSolver::Run()
{
  GlobalData::GetInstance()->m_cycle += 1;
  GlobalData::GetInstance()->m_time += m_dTime;

  // Get the Displacement, Velocity, Acceleration
  Vec &disp = GlobalData::GetInstance()->m_state;
  Vec &velo = GlobalData::GetInstance()->m_velo;
  Vec &acce = GlobalData::GetInstance()->m_acce;

  // Get the Internal Force and External Force
  Vec &fint = GlobalData::GetInstance()->m_fint;
  Vec &fhat = GlobalData::GetInstance()->m_fhat;

  VecAXPY(velo, 0.5*m_dTime, acce);
  VecAxPy(disp, m_dTime, velo);
  
  GlobalData::GetInstance()->m_elements->AssembleInternalForce(fint);

  double lam = 0;
  GlobalData::GetInstance()->m_dofs->SetConstrainFactor(lam);

  Vec res;
  VecDuplicate(fhat, &res);
  VecCopy(fhat, res);
  VecScale(res, lam);
  VecAXPY(res, -1., fint);

  KSP ksp,
  KSPCreate(ksp, PETSC_COMM_WORLD);

  GlobalData::GetInstance()->m_dofs->Solve(m_lumpedMass, res, acce, ksp);

  VecAXPY(velo, 0.5*m_dTime, acce);

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
  MatMult(m_lumpedMass, velo, tempVec);
  VecDot(tempVec,velo, &energy);
  
  std::cout << std::setw(7) << GlobalData::GetInstance()->m_cycle;

  std::cout.precision(3);
  std::cout << std::setw(13) << std::setiosflags(std::ios::scientific)
            << GlobalData::GetInstance()->m_time
            << std::setw(13) << std::setiosflags(std::ios::scientific) << energy << std::endl;
}