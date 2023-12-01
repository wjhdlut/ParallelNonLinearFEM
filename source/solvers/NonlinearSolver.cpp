/**
 * @File Name:     NonlinearSolver.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <iostream>
#include <regex>
#include <fstream>
#include <solvers/NonlinearSolver.h>
#include <util/DataStructure.h>
#include <petscksp.h>

NonlinearSolver::NonlinearSolver(const nlohmann::json &props)
                 : BaseModule(props)
{
  Initialize(props);
}

NonlinearSolver::~NonlinearSolver()
{}

void NonlinearSolver::Initialize(const nlohmann::json &props)
{
  GlobalData::GetInstance()->m_lam = 0.;

  // Two Type Of Method to Read Solver Information
  std::string fileName = props.at("input");
  ReadData(fileName);

  if(0 == m_maxCycle && m_myProps.contains("maxCycle")){
    Tools::GetParameter(m_maxCycle, "maxCycle", m_myProps);
    m_iterMaxVec.assign(m_maxCycle, 10);
    m_tolVec.assign(m_maxCycle, 1.0e-3);
  }

  if(0 == m_maxCycle)
    throw "Please Check Solver Information and Assign the Number Of Load Step";
}

void NonlinearSolver::Run()
{
  int &numOfCycle = GlobalData::GetInstance()->m_cycle;
  GlobalData::GetInstance()->m_cycle += 1;
  
  // Update the Load Factor
  UpdateLoadFactor(numOfCycle);

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
  while(error > m_tolVec[numOfCycle - 1])
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
    
    if(GlobalData::GetInstance()->m_iiter == m_iterMaxVec[numOfCycle - 1])
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

void NonlinearSolver::ReadData(const std::string &fileName)
{
  std::ifstream fin(fileName, std::ios::in);
  if(!fin.is_open()){
    std::cout << fileName << " open failed!!" << std::endl;
    exit(-1);
  }

  std::string line = "";
  std::regex pattern("[\\s]{2,}");

  while(true)
  {
    getline(fin, line);
    if(line.npos != line.find("\r")) line.erase(line.find("\r"));

    if(line.npos != line.find("<NonlinearSolver>"))
    {
      while(true)
      {
        getline(fin, line);
        if(line.npos != line.find("\r"))
          line.erase(line.find("\r"));
        
        if(line.npos != line.find("</NonlinearSolver>")){
           m_maxCycle = m_tolVec.size();
           fin.close();
          return;
        }

        if(0 == line.size()) continue;

        line = std::regex_replace(line, pattern, " ");
        std::vector<std::string> strData = Tools::StringSplit(line, ";");
        for(auto iter = strData.begin(); iter != strData.end() - 1; iter++)
        {
          std::string tempStr = Tools::StringStrip(*iter);
          std::vector<std::string> b = Tools::StringSplit(tempStr, " ");
          
          // Skip the Comment Line
          if("//" == b[0].substr(0, 2) || "#" == b[0].substr(0, 1)) break;
          
          // Read Incremental Load Factor
          m_lamVec.emplace_back(std::stod(b[0]));
          
          // Read Convergence Tolerence
          m_tolVec.emplace_back(std::stod(b[1]));
          
          // Read Max. No. of Iterations
          m_iterMaxVec.emplace_back(std::stoi(b[2]));

          // Output Control Flags for Results TO DO WORK
        }
      }
    }
    if(fin.eof()) break;
  }
  fin.close(); 
}

void NonlinearSolver::UpdateLoadFactor(const int &numOfCycle)
{
  if(0 == m_lamVec.size()){
    GlobalData::GetInstance()->m_lam = 1.0 * GlobalData::GetInstance()->m_cycle;
  }
  else{
    GlobalData::GetInstance()->m_lam += m_lamVec[numOfCycle - 1];
  }
}