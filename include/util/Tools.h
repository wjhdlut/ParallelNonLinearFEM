/**
 * @File Name:     Tools.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <petscvec.h>

namespace Tools
{
  /**
   * @Brief: Split a string into a list where each word is a list item
   * 
   * @param strData 
   * @param tag 
   * @param maxSplit 
   * @return std::vector<std::string> 
   */
  std::vector<std::string> StringSplit(std::string &strData,
                                       const std::string &tag,
                                       const int maxSplit = 20);

  /**
   * @Brief: Remove the Assignment CHaracter at the beginning and at the end of the string
   * 
   * @param strData 
   * @param tag 
   * @return std::string 
   */
  std::string StringStrip(const std::string &strData, const std::string &tag = " ");
  
  /**
   * @Brief: Get the Variable Type
   * 
   * @param value 
   * @return std::string 
   */
  std::string GetVarType(const std::string&value);

  /**
   * @Brief:  Get the Int Value with Name in Props JSON Object
   * 
   * @param value 
   * @param name 
   * @param props 
   */
  void GetParameter(int &value, const std::string &name, const nlohmann::json &props);

  /**
   * @Brief: Get the Double Value with Name in Props JSON Object
   * 
   * @param value 
   * @param name 
   * @param props 
   */
  void GetParameter(double &value, const std::string &name, const nlohmann::json &props);

  /**
   * @Brief: Get the Bool Value with Name in Props JSON Object
   * 
   * @param value 
   * @param name 
   * @param props 
   */
  void GetParameter(bool &value, const std::string &name, const nlohmann::json &props);

  /**
   * @Brief: Get the String Value with Name in Props JSON Object
   * 
   * @param value 
   * @param name 
   * @param props 
   */
  void GetParameter(std::string &value, const std::string &name, const nlohmann::json &props);

  /**
   * @Brief:  Print PETSc Vec into File
   * 
   * @param data 
   * @param fileName 
   * @return PetscErrorCode 
   */
  PetscErrorCode PrintVecIntoFile(const Vec&data, const std::string &fileName);
}

#endif // TOOLS_H