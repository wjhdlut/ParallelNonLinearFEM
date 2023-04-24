#include <solvers/Solver.h>
#include <util/DataStructure.h>
#include <util/ObjectFactory.h>
// #include <iostream>

Solver::Solver()
{
  nlohmann::json &solverProps = GlobalData::GetInstance()->m_props.at("solver");

  const std::string&solveType = solverProps.at("type");

  GlobalData::GetInstance()->m_props["currentModule"] = "solver";

  m_solver = ObjectFactory::CreateObject<BaseModule>(solveType, GlobalData::GetInstance()->m_props);
}

Solver::~Solver()
{
}

void Solver::Run()
{
  if(nullptr != m_solver) m_solver->Run();
}
