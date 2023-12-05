#include <fstream>
#include <iostream>
#include <variant>
#include <regex>

#include "../include/FileParse.h"
#include "../include/DataStructure.h"

void StoreValue(nlohmann::json&db, const std::string&key, const std::vector<std::string> &value)
{
  if(!db.contains(key)) db[key] = nlohmann::json::array();
  for(auto iValue : value)
  {
    db[key].emplace_back(iValue);
  }
}


void StoreValue(nlohmann::json&db, const std::string&key, const std::string &value)
{
  std::regex intPattern("[\\d]+");
	std::regex floatPattern("[\\d]+[.][\\d]*[eE]?[+-]?[\\d]*");
	std::smatch result;

  if("true" == value) db[key] = true;
	else if("false" == value) db[key] = false;
  else if(regex_match(value, result, intPattern)) db[key] = std::stoi(value);
	else if(regex_match(value, result, floatPattern)) db[key] = std::stod(value);
  else db[key] = value;
}

std::string ReadItem(std::vector<std::string>&strVec, nlohmann::json&db)
{
  std::vector<std::string> l2;
  if(strVec[0].npos != strVec[0].find(".")){
    l2 = Tools::StringSplit(strVec[0], ".", 1);
    return 0;
  }
  else{
    l2 = Tools::StringSplit(strVec[1], ";", 1);
    
    if("[" == l2[0].substr(0, 1)){
      int tempInt = std::string::npos != l2[0].find("]") ? 2 : 1;
      std::string tempStr = l2[0].substr(1, l2[0].size() - tempInt);
      std::vector<std::string> l3 = Tools::StringSplit(tempStr, ",");
      StoreValue(db, strVec[0], l3);

      std::vector<std::string> l4 = Tools::StringSplit(l2[1], "]", 1);
      if (l4.size() == 2){
        if(";" == l4[1].substr(0,1)) l4[1].erase(0, 1);
        l2 = l4;
        while(2 == l3.size() && 0 != l3[1].size()){
          l3 = Tools::StringSplit(l2[0], ";", 1);
          l4 = Tools::StringSplit(l3[0], ",", 1);
          StoreValue(db, strVec[0], l4);
        }
      }
    }
    else
      StoreValue(db, strVec[0], l2[0]);
    return l2[1];
  }
}

std::string ReadBlock(nlohmann::json &db, const std::string &ln)
{
  std::string fileData = ln;
  while(true)
  {
    if("include" == fileData.substr(0, 7))
    {
      std::vector<std::string> l1 = Tools::StringSplit(fileData, ";", 1);
    }

    std::vector<std::string> l1 = Tools::StringSplit(fileData, "=", 1);
    
    if(1 == l1.size()) return fileData;

    if("};" == l1[0].substr(0, 2)) return l1[0].erase(0,2) + "=" + l1[1];

    if("{" == fileData.substr(0, 1)){
      nlohmann::json child = nlohmann::json::object();
      fileData.erase(0, 1);
      fileData = ReadBlock(child, fileData);
      db[l1[0]] = child;
    }
    else{
      fileData = ReadItem(l1, db);
    }
  }
}

void FileParse(nlohmann::json &db, const std::string &fileName)
{
  std::ifstream inputFile(fileName, std::ios::in);
  if(!inputFile.is_open()){
    std::cout << fileName << " open failed!!" << std::endl;
    exit(-1);
  }
  
  std::string ln = "";
  std::string line = "";
  while(!inputFile.eof())
  {
    getline(inputFile, line);
    
    if(line.size() > 0 && line.substr(line.find_first_not_of(" "), 1) != "#")
    {
      if(line.npos != line.find_last_of("\r"))
        line.erase(line.find_last_of("\r"));
      ln.append(line);
    }
  }
  while (ln.npos != ln.find(" ")) ln.erase(ln.find(" "), 1);
  while (ln.npos != ln.find("\"")) ln.erase(ln.find("\""), 1);

  ReadBlock(db, ln);
  inputFile.close();

  // std::cout << std::setw(4) << db << std::endl;
}