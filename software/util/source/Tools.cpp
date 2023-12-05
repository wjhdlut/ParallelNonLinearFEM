/**
 * @File Name:     Tools.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <regex>
#include <iostream>

#include "../include/Tools.h"

namespace Tools
{
std::vector<std::string> StringSplit(std::string &strData, const std::string &tag, const int maxSplit)
{
  int count = 0;
  std::vector<std::string> strVec;
  do
  {
    if(strData.npos == strData.find_first_of(tag)) break;
    std::string l1 = strData.substr(0, strData.find_first_of(tag));
    strData.erase(0, strData.find_first_of(tag) + 1);
    // if(strData.size() > 0) strData = StringStrip(strData);
    strVec.emplace_back(l1);
    count += 1;
  }while(count < maxSplit);
  strVec.emplace_back(strData);
  return strVec;
}

std::string StringStrip(const std::string &strData, const std::string &tag)
{
  int tagLength = tag.size();
  std::string str = strData;
  if(tag == str.substr(0, tagLength)) str.erase(0, tagLength);
  
  int strDataLength = str.size();
  if(tag == str.substr(strDataLength-tagLength, tagLength)) str.erase(strDataLength-tagLength);

  return str;
}

std::string GetVarType(const std::string&value)
{
  std::regex intPattern("[\\d]+");
	std::regex floatPattern("[\\d]+[.][\\d]*");
  std::smatch result;

  if("true" == value || "false" == value) return "bool";
  else if(regex_match(value, result, intPattern)) return "int";
	else if(regex_match(value, result, floatPattern)) return "double";
  else return "string";
}

void GetParameter(int &value, const std::string &name, const nlohmann::json &props)
{
  if(props.contains(name))
  {
    if(props.at(name).is_string())
    {
      std::string temp = props.at(name);
      value = std::stoi(temp);
    }
    else
    {
      value = props.at(name);
    }
  }
}

void GetParameter(double &value, const std::string &name, const nlohmann::json &props)
{
  if(props.contains(name))
  {
    if(props.at(name).is_string())
    {
      std::string temp = props.at(name);
      value = std::stod(temp);
    }
    else
    {
      value = props.at(name);
    }
  }
}

void GetParameter(bool &value, const std::string &name, const nlohmann::json &props)
{
  if(props.contains(name))
    value = props.at(name);
}

void GetParameter(std::string &value, const std::string &name, const nlohmann::json &props)
{
  if(props.contains(name))
    value = props.at(name);
}

PetscErrorCode PrintVecIntoFile(const Vec&data, const std::string&fileName)
{
  PetscErrorCode ierr;
  FILE *fr;
  fr=fopen(fileName.c_str(), "w");
  
  int numOfVecSize;
  double temp;
  ierr = VecGetSize(data, &numOfVecSize);
  for(int index = 0; index < numOfVecSize; index++){
    ierr = VecGetValues(data, 1, &index, &temp);
    ierr = PetscFPrintf(PETSC_COMM_WORLD, fr, "%1.6e\n", temp);
  }
  fclose(fr);

  return ierr;
}
}