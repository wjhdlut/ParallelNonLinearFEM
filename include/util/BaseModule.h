#ifndef BASEMOULE_H
#define BASEMOULE_H

#include <nlohmann/json.hpp>
#include <iostream>
#include <util/Tools.h>

class BaseModule
{
public:
  inline BaseModule(const nlohmann::json &props){
    std::string currentModule = "solver";
    if(props.contains("currentModule")){
      if(props.contains(props.at("currentModule"))){
        currentModule = props.at("currentModule");
      }
    }

    if(props.contains(currentModule)) m_myProps = props.at(currentModule);
  }
  virtual ~BaseModule() = default;

  virtual void Run() = 0;

protected:
  inline void GetParameter(int &value, const std::string &name){
    Tools::GetParameter(value, name, m_myProps);
  }

  inline void GetParameter(double &value, const std::string &name){
    Tools::GetParameter(value, name, m_myProps);
  }

  inline void GetParameter(std::string &value, const std::string &name){
    Tools::GetParameter(value, name, m_myProps);
  }

  inline void GetParameter(bool &value, const std::string &name){
    Tools::GetParameter(value, name, m_myProps);
  }

protected:
  nlohmann::json m_myProps;
};

#endif // BASEMOULE_H