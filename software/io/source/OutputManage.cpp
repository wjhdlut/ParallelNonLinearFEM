/**
 * @File Name:     OutputManage.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <string>

#include "../include/OutputManager.h"
#include "../../util/include/DataStructure.h"
#include "../../util/include/ObjectFactory.h"

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

    std::shared_ptr<BaseModule> ins = ObjectFactory::CreateObject<BaseModule>(ioType, GlobalData::GetInstance()->m_props);
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