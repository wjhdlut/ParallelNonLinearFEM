/**
 * @File Name:     OutputWriter.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/OutputWriter.h"
#include "../../util/include/DataStructure.h"

OutputWriter::OutputWriter(const nlohmann::json &props) : BaseModule(props)
{
  
}

OutputWriter::~OutputWriter()
{}

void OutputWriter::Run()
{
  GlobalData::GetInstance()->PrintNodes();
}