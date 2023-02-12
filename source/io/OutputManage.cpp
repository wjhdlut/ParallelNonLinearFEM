
#include <string>
#include <util/DataStructure.h>
#include <io/OutputManager.h>
#include <util/ObjectFactory.h>

OutputManager::OutputManager()
{
  std::vector<std::string> outputModules = GlobalData::GetInstance()->m_props["outputModules"];
  for (auto name : outputModules)
  {
    GlobalData::GetInstance()->m_props["currentModule"] = name;
    std::string ioType = name;
    if (GlobalData::GetInstance()->m_props.contains(name))
    {
      nlohmann::json &moduleProps = GlobalData::GetInstance()->m_props.at(name);
      if (moduleProps.contains("type"))
        ioType = moduleProps.at("type");
    }

    std::shared_ptr<BaseModule> ins = ObjectFactory::CreateObject<BaseModule>(ioType);
    if (nullptr != ins) m_outman.emplace_back(ins);
  }
}

OutputManager::~OutputManager()
{}

void OutputManager::Run()
{
  for(auto output : m_outman)
    output->Run();
}