/**
 * @File Name:     GraphWriter.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <io/GraphWriter.h>
#include <util/DataStructure.h>
#include <iostream>

GraphWriter::GraphWriter(const nlohmann::json &props) : BaseModule(props)
{ 
  std::string fileName = GlobalData::GetInstance()->m_prefix + ".out";
  std::cout << GlobalData::GetInstance()->m_prefix << std::endl;
  
  m_outFileStream.open(fileName, std::ios::out);
  m_columns = m_myProps.at("columns");

  for(auto iColumns : m_columns)
    m_outFileStream << iColumns << "       ";
  m_outFileStream << std::endl;
}

GraphWriter::~GraphWriter()
{
}

void GraphWriter::Run()
{
  VectorXd data = VectorXd::Zero(m_columns.size());
  for(int i = 0; i < m_columns.size(); i++){
    const nlohmann::json &columnsData = m_myProps.at(m_columns[i]);
    std::string type;
    Tools::GetParameter(type, "type", columnsData);

    int nodeID;
    Tools::GetParameter(nodeID, "node", columnsData);
    if(GlobalData::GetInstance()->m_outputName.end() != 
       std::find(GlobalData::GetInstance()->m_outputName.begin(),
       GlobalData::GetInstance()->m_outputName.end(), type))  
    {
      data = GlobalData::GetInstance()->GetData(type, nodeID);
    }
    else{
      Vec tempV = GetGlobalData(type);
      dofType = columnsData.at("dof");
      index = GlobalData::GetInstance()->m_dofs->GetForType(nodeID, dofType);
      VecGetValues(tempV, 1, &index, &tempDouble);
      Tools::GetParameter(factor, "factor", columnsData);
      data(i) = tempDouble * factor;
    }
    m_outFileStream << data(i) << "       ";
  }
  m_outFileStream << std::endl;

  m_data.emplace_back(data);
}

Vec& GraphWriter::GetGlobalData(const std::string &name)
{
  if("state" == name){
    return (GlobalData::GetInstance()->m_state);
  }
  else if("fint" == name){
    return (GlobalData::GetInstance()->m_fint);
  }
  else{
    throw "no output data with name: " + name;
  }
}