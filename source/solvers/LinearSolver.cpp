/**
 * @File Name:     LinearSolver.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <solvers/LinearSolver.h>
#include <petscksp.h>
#include <util/DataStructure.h>
#include <iostream>

LinearSolver::LinearSolver(const nlohmann::json&props) : BaseModule(props)
{}

LinearSolver::~LinearSolver()
{}

void LinearSolver::Run()
{
  GlobalData::GetInstance()->m_cycle += 1;

  Vec &a = GlobalData::GetInstance()->m_state;
  Vec &fext = GlobalData::GetInstance()->m_fhat;
  Vec &fint = GlobalData::GetInstance()->m_fint;

  Mat K;
  GlobalData::GetInstance()->m_elements->AssembleTangentStiffness(K, fint);
  // std::cout << "Stiffness Matrix = " << std::endl;
  // MatView(K, PETSC_VIEWER_STDOUT_WORLD);

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  GlobalData::GetInstance()->m_dofs->Solve(K, fext, a, ksp);
  // std::cout << "displacement = " << std::endl;
  // VecView(a, PETSC_VIEWER_STDOUT_WORLD);

  GlobalData::GetInstance()->m_elements->AssembleInternalForce(fint);

  GlobalData::GetInstance()->m_elements->CommitHistory();
  
  GlobalData::GetInstance()->m_active = false;

  MatDestroy(&K);
  KSPDestroy(&ksp);
}