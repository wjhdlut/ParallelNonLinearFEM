#ifndef FILEPARSE_H
#define FILEPARSE_H

#include <util/Tools.h>
#include <nlohmann/json.hpp>


/**
 * @Brief:   store related data into json
 * 
 * @param db                    [in/out]    json data
 * @param key                   [in]        key in json data
 * @param value                 [in]        value in json data
 */
void StoreValue(nlohmann::json&db, const std::string&key, const std::vector<std::string> &value);

/**
 * @Brief:  store relate data into json
 * 
 * @param db                    [in/out]    json data
 * @param key                   [in]        key in json data
 * @param value                 [in]        value in json data
 */
void StoreValue(nlohmann::json&db, const std::string&key, const std::string &value);

/**
 * @Brief:  read data in input file with []
 * 
 * @param strData 
 * @param bd 
 * @return std::string 
 */
std::string ReadItem(std::vector<std::string> &strData, nlohmann::json &bd);


/**
 * @Brief:  read data in input file with {}
 * 
 * @param db 
 * @param ln 
 * @return std::string 
 */
std::string ReadBlock(nlohmann::json &db, const std::string &ln);

/**
 * @Brief: Parse input file
 * 
 * @param db 
 * @param fileName 
 */
void FileParse(nlohmann::json &db, const std::string&fileName);

#endif // FILEPARSE_H