/**
 * @File Name:     RiksSolver.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <solvers/RiksSolver.h>
#include <util/DataStructure.h>
#include <iostream>

RiksSolver::RiksSolver(const nlohmann::json &props) : BaseModule(props)
{
  Initialize(props);
}

RiksSolver::~RiksSolver()
{
}

void RiksSolver::Initialize(const nlohmann::json &props)
{
  GlobalData::GetInstance()->m_lam = 1.0;
  
  Tools::GetParameter(m_fixedStep, "fixedStep", m_myProps);
  Tools::GetParameter(m_maxLam, "maxLam", m_myProps);
}

void RiksSolver::Run()
{
  GlobalData::GetInstance()->m_cycle += 1;

  Vec &a = GlobalData::GetInstance()->m_state;
  Vec &Da = GlobalData::GetInstance()->m_Dstate;
  Vec &fhat = GlobalData::GetInstance()->m_fhat;
  Vec &fint = GlobalData::GetInstance()->m_fint;

  Vec fext, da1, res, daPrev;
  VecDuplicate(fhat, &fext);
  VecCopy(fhat, fext);
  VecDuplicate(fhat, &daPrev);
  VecCopy(Da, daPrev);
  // std::cout << "fext: " << std::endl;
  // VecView(fext, PETSC_VIEWER_STDOUT_WORLD);

  VecDuplicate(fhat, &da1);
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
  if(1 == GlobalData::GetInstance()->m_cycle)
  {
    GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);
    // std::cout << "StiffnessMatrix:" << std::endl;
    // MatView(K, PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "fext = " << std::endl;
    // VecView(fext, PETSC_VIEWER_STDOUT_WORLD);
    VecScale(fext, GlobalData::GetInstance()->m_lam);
    GlobalData::GetInstance()->m_dofs->Solve(K, fext, da1, ksp);
    // std::cout << "da1: " << std::endl;
    // VecView(da1, PETSC_VIEWER_STDOUT_WORLD);
    dLam = GlobalData::GetInstance()->m_lam;
  }
  else
  {
    VecCopy(daPrev, da1);
    // std::cout << "da1: " << std::endl;
    // VecView(da1, PETSC_VIEWER_STDOUT_WORLD);
    VecScale(da1, m_factor);

    dLam = m_factor * m_dLamPrev;
    GlobalData::GetInstance()->m_lam += dLam;
  }

  VecAXPY(a, 1., da1);
  VecCopy(da1, Da);
  
  GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);
  // std::cout << "Stiffness Matrix:" << std::endl;
  // MatView(K, PETSC_VIEWER_STDOUT_WORLD);
  // std::cout << "fint: " << std::endl;
  // VecView(fint, PETSC_VIEWER_STDOUT_WORLD);

  VecScale(fext, GlobalData::GetInstance()->m_lam);
  VecWAXPY(res, -1., fint, fext);

  Vec tempVecD1, tempVecD2, dDa;
  VecDuplicate(fext, &tempVecD1);
  VecDuplicate(fext, &tempVecD2);
  VecDuplicate(fext, &dDa);

  double dDLam = 0., tempDValue1= 0., tempDValue2 = 0.;

  // Loop For Iterative Step
  while(error > m_tol)
  {
    GlobalData::GetInstance()->m_iiter += 1;
    
    GlobalData::GetInstance()->m_dofs->Solve(K, fhat, tempVecD1, ksp);
    // std::cout << "fext:" << std::endl;
    // VecView(fext, PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "tempVecD1: " << std::endl;
    // VecView(tempVecD1, PETSC_VIEWER_STDOUT_WORLD);
    GlobalData::GetInstance()->m_dofs->Solve(K, res, tempVecD2, ksp);
    // std::cout << "res:" << std::endl;
    // VecView(res, PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "tempVecD2:" << std::endl;
    // VecView(tempVecD2, PETSC_VIEWER_STDOUT_WORLD);

    VecDot(da1, tempVecD2, &tempDValue2);
    VecDot(da1, tempVecD1, &tempDValue1);
    dDLam = - tempDValue2 / tempDValue1;

    VecWAXPY(dDa, dDLam, tempVecD1, tempVecD2);

    dLam += dDLam;
    GlobalData::GetInstance()->m_lam += dDLam;

    VecAXPY(a, 1., dDa);
    VecAXPY(Da, 1., dDa);
    // std::cout << "a=" << std::endl;
    // VecView(a, PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "Da=" << std::endl;
    // VecView(Da, PETSC_VIEWER_STDOUT_WORLD);

    GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);
    // std::cout << "StiffnessMatrix:" << std::endl;
    // MatView(K, PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "fint: " << std::endl;
    // VecView(fint, PETSC_VIEWER_STDOUT_WORLD);

    VecCopy(fhat, fext);
    VecScale(fext, GlobalData::GetInstance()->m_lam);
    VecWAXPY(res, -1., fint, fext);

    // std::cout << "res: " << std::endl;
    // VecView(res, PETSC_VIEWER_STDOUT_WORLD);

    GlobalData::GetInstance()->m_dofs->Norm(res, tempDValue1);
    GlobalData::GetInstance()->m_dofs->Norm(fext, tempDValue2);
    error = tempDValue1 / tempDValue2;

    // Print Iterative Information
    std::cout.precision(2);
    std::cout << "  Iter " << GlobalData::GetInstance()->m_iiter
              << " : error = " << std::setiosflags(std::ios::scientific) << error << std::endl;

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

  m_dLamPrev = dLam;

  if(GlobalData::GetInstance()->m_lam > m_maxLam ||
     GlobalData::GetInstance()->m_cycle > 1000)
     GlobalData::GetInstance()->m_active = false;

  // free the space
  VecDestroy(&fext);
  VecDestroy(&da1);
  VecDestroy(&res);
  VecDestroy(&tempVecD1);
  VecDestroy(&tempVecD2);
  VecDestroy(&dDa);
  VecDestroy(&daPrev);
  MatDestroy(&K);
  KSPDestroy(&ksp);
}