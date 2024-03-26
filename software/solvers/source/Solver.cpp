/**
 * @File Name:     Solver.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/Solver.h"
#include "../../util/include/DataStructure.h"
#include "../../util/include/ObjectFactory.h"

Solver::Solver()
{
  Initialize();
}

Solver::~Solver()
{
}

void Solver::Initialize()
{
  nlohmann::json &solverProps = GlobalData::GetInstance()->m_props.at("solver");

  const std::string&solveType = solverProps.at("type");

  GlobalData::GetInstance()->m_props["currentModule"] = "solver";

  m_solver = ObjectFactory::CreateObject<BaseModule>(solveType, GlobalData::GetInstance()->m_props);
  if(nullptr == m_solver)
  {
    std::cout << solveType << " Solver Created Failed!!!" << std::endl;
    exit(-1);
  }

  GlobalData::GetInstance()->m_dofs->CompBasicVariable();
}

void Solver::Run()
{
  if(nullptr != m_solver) m_solver->Run();
}
