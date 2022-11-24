
#include <util/Tools.h>
#include <regex>

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
}