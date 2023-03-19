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
  std::vector<double> data;
  for(auto iColumns : m_columns){
    const nlohmann::json &columnsData = m_myProps.at(iColumns);
    std::string type;
    Tools::GetParameter(type, "type", columnsData);

    int nodeID;
    Tools::GetParameter(nodeID, "node", columnsData);
    if(GlobalData::GetInstance()->m_outputName.end() != 
       std::find(GlobalData::GetInstance()->m_outputName.begin(),
       GlobalData::GetInstance()->m_outputName.end(), type)){
      data = GlobalData::GetInstance()->GetData(type, nodeID);
    }
    else{
      Vec tempV = GetGlobalData(type);
      dofType = columnsData.at("dof");
      index = GlobalData::GetInstance()->m_dofs->GetForType(nodeID, dofType);
      VecGetValues(tempV, 1, &index, &tempDouble);
      Tools::GetParameter(factor, "factor", columnsData);
      tempDouble *= factor;
      data.emplace_back(tempDouble);
    }
    m_outFileStream << tempDouble << "       ";
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