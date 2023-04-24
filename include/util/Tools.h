#ifndef TOOLS_H
#define TOOLS_H


#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <petscvec.h>

namespace Tools
{
  std::vector<std::string> StringSplit(std::string &strData, const std::string &tag, const int maxSplit = 10);

  std::string StringStrip(const std::string &strData, const std::string &tag = " ");

  std::string GetVarType(const std::string&value);

  void GetParameter(int &value, const std::string &name, const nlohmann::json &props);

  void GetParameter(double &value, const std::string &name, const nlohmann::json &props);

  void GetParameter(bool &value, const std::string &name, const nlohmann::json &props);

  void GetParameter(std::string &value, const std::string &name, const nlohmann::json &props);

  PetscErrorCode PrintVecIntoFile(const Vec&data, const std::string &fileName);
}

#endif // TOOLS_H