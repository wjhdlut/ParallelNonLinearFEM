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
    ierr = PetscFPrintf(PETSC_COMM_WORLD, fr, "%1.20e\n", temp);
  }
  fclose(fr);

  return ierr;
}

PetscErrorCode PrintMatIntoFile(const Mat &data, const std::string &fileName)
{
  PetscErrorCode ierr;
  FILE *fr;
  fr=fopen(fileName.c_str(), "w");
  
  int rowsOfMat, colsOfMat;
  double temp;
  ierr = MatGetSize(data, &rowsOfMat, &colsOfMat);
  for(int row = 0; row < rowsOfMat; row++){
    for(int col = 0; col < colsOfMat; col++){
      ierr = MatGetValue(data, row, col, &temp);
      ierr = PetscFPrintf(PETSC_COMM_WORLD, fr, "%1.20e   ", temp);
    }
    ierr = PetscFPrintf(PETSC_COMM_WORLD, fr, "\n");
  }
  fclose(fr);

  return ierr;
}
void tempFunction(const Mat &K)
{
  Mat tempTransMat;
  MatCreate(PETSC_COMM_WORLD, &tempTransMat);
  MatSetSizes(tempTransMat, PETSC_DECIDE, PETSC_DECIDE, 90, 90);
  MatSetFromOptions(tempTransMat);
  MatSetUp(tempTransMat);
  std::vector<int> temp = { 1,  2, 37, 38, 21, 22, 89, 90,  7,  8, 57, 58,  5,  6, 65, 66,  9, 10, 75, 76, 15, 16, 43, 44, 39, 40, 55, 56, 51, 52, 67, 68, 77, 78, 25, 26, 41, 42,
                           31, 32, 85, 86, 33, 34, 53, 54, 35, 36, 69, 70, 27, 28, 79, 80, 29, 30, 49, 50, 45, 46, 63, 64, 59, 60, 71, 72, 81, 82,  3,  4, 47, 48, 19, 20,
                           87, 88, 23, 24, 61, 62, 13, 14, 73, 74, 11, 12, 83, 84, 17, 18};
  for(int i = 0; i < 90; i++)
  {
    MatSetValue(tempTransMat, i, temp[i]-1, 1., INSERT_VALUES);
  }
  MatAssemblyBegin(tempTransMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(tempTransMat, MAT_FINAL_ASSEMBLY);
  
  Mat tempK;
  MatTransposeMatMult(tempTransMat, K, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tempK);
  // MatMatMult(tempK, tempTransMat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &K);
  std::cout << "stiffness Matrix 2 = " << std::endl;
  MatView(K, PETSC_VIEWER_STDOUT_WORLD);
}
}