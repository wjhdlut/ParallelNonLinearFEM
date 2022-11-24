#ifndef FILEPARSE_H
#define FILEPARSE_H

#include <util/Tools.h>
#include <nlohmann/json.hpp>

void StoreValue(nlohmann::json&db, const std::string&key, const std::vector<std::string> &value);

void StoreValue(nlohmann::json&db, const std::string&key, const std::string &value);

std::string ReadItem(std::vector<std::string> &strData, nlohmann::json &bd);

std::string ReadBlock(nlohmann::json &db, const std::string &ln);

void FileParse(nlohmann::json &db, const std::string&fileName);

#endif // FILEPARSE_H