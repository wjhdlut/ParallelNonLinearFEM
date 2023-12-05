/**
 * @File Name:     BaseModule.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef BASEMOULE_H
#define BASEMOULE_H

#include <iostream>

#include "Tools.h"
#include "../../nlohmann/json.hpp"

class BaseModule
{
public:
  /**
   * @Brief: Construct a new Base Module object
   * 
   * @param props 
   */
  inline BaseModule(const nlohmann::json &props){
    std::string currentModule = "solver";
    if(props.contains("currentModule")){
      if(props.contains(props.at("currentModule"))){
        currentModule = props.at("currentModule");
      }
    }

    if(props.contains(currentModule)) m_myProps = props.at(currentModule);
  }

  /**
   * @Brief: Destroy the Base Module object
   * 
   */
  virtual ~BaseModule() = default;

  /**
   * @Brief: Run the Instance
   * 
   */
  virtual void Run() = 0;

  /**
   * @Brief: Read Data From File
   * 
   * @param fileName 
   */
  inline virtual void ReadData(const std::string &fileName) {}

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

  inline virtual void Initialize(const nlohmann::json &props) {

  }

protected:
  nlohmann::json m_myProps;
};

#endif // BASEMOULE_H