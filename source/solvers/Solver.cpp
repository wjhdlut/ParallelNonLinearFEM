#include <solvers/Solvers.h>
#include <util/DataStructure.h>
#include <util/ObjectFactory.h>

Solvers::Solvers()
{
  nlohmann::json &solverProps = GlobalData::GetInstance()->m_props.at("solver");

  const std::string&solveType = solverProps.at("type");

  GlobalData::GetInstance()->m_props["currentModule"] = "solver";

  m_solver = ObjectFactory::CreateObject<BaseModule>(solveType, GlobalData::GetInstance()->m_props);
}

Solvers::~Solvers()
{}

void Solvers::Run()
{
  m_solver->Run();
}
