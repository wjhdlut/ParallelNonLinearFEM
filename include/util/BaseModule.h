#ifndef BASEMOULE_H
#define BASEMOULE_H

#include <nlohmann/json.hpp>
#include <iostream>

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
    if(m_myProps.contains(name))
      if(m_myProps.at(name).is_string()){
        std::string temp = m_myProps.at(name);
        value = std::stoi(temp);
      }
      else
        value = m_myProps.at(name);
    else
      throw "no properties with name: " + name;
  }

  inline void GetParameter(double &value, const std::string &name){
    if(m_myProps.contains(name))
      if(m_myProps.at(name).is_string()){
        std::string temp = m_myProps.at(name);
        value = std::stod(temp);
      }
      else
        value = m_myProps.at(name);
    else
      throw "no properties with name: " + name;
  }

  inline void GetParameter(std::string &value, const std::string &name)
  {
    if(m_myProps.contains(name))
      value = m_myProps.at(name);
    else
      throw "no properties with name: " + name;
  }

  inline void GetParameter(bool &value, const std::string &name)
  {
    if(m_myProps.contains(name))
      value = m_myProps.at(name);
    else
      throw "no properties with name: " + name;
  }

protected:
  nlohmann::json m_myProps;
};

#endif // BASEMOULE_H