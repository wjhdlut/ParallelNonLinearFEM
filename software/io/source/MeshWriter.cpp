/**
 * @File Name:     MeshWriter.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <iostream>
#include <fstream>

#include "../include/MeshWriter.h"
#include "../../util/include/ShapeFunctions.h"
#include "../../util/include/DataStructure.h"

MeshWriter::MeshWriter(const nlohmann::json &props) : BaseModule(props)
{
  GetParameter(m_interval, "interval");
}

MeshWriter::~MeshWriter()
{}

void MeshWriter::Run()
{
  if(0 != (GlobalData::GetInstance()->m_cycle%m_interval)) return;
  
  // output vtk file
  WriteVTKFile();

  // write pvd file
  WritePVDFile();

  m_k += 1;
}

void MeshWriter::WritePointInfo(std::fstream &vtkfile)
{
  vtkfile << "<PointData>" << std::endl;
  WritePointDisplaceData(vtkfile);

  // output stress information
  WritePointStressData(vtkfile);

  vtkfile << "</PointData>" << std::endl;
  vtkfile << "<CellData>" << std::endl;
  vtkfile << "</CellData>" << std::endl;

  // output node coordinate information
  WritePointCoordsData(vtkfile);
}

void MeshWriter::WritePointDisplaceData(std::fstream &vtkfile)
{
  vtkfile << "<DataArray type=\"Float64\" Name=\"displacement\""
          << " NumberOfComponents=\"3\" format=\"ascii\" >" << std::endl;

  // output node displacement information
  double nodeCoordValue = 0.;
  int dofIndex = 0;
  for(auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords){
    for(auto dofType : GlobalData::GetInstance()->m_dofs->GetDofType()){
      dofIndex = GlobalData::GetInstance()->m_dofs->GetForType(node.first, dofType);
      VecGetValues(GlobalData::GetInstance()->m_state, 1, &dofIndex, &nodeCoordValue);
      vtkfile << nodeCoordValue << "  ";
    }
    if(node.second.size() < 3) vtkfile << "  0.";
    vtkfile << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;
}

void MeshWriter::WritePointStressData(std::fstream &vtkfile)
{
  // output sigma_xx
  MatrixXd stress = GlobalData::GetInstance()->GetData("stresses");
  vtkfile << "<DataArray type=\"Float64\" Name=\"sigma_xx\""
          << " NumberOfComponents=\"1\" format=\"ascii\" >" << std::endl;
  for (auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords)
  {
    int index = GlobalData::GetInstance()->m_dofs->GetIndex(node.first);
    vtkfile << stress(index, 0) << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;

  // output sigma_yy
  vtkfile << "<DataArray type=\"Float64\" Name=\"sigma_yy\""
          << " NumberOfComponents=\"1\" format=\"ascii\" >" << std::endl;
  for (auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords)
  {
    int index = GlobalData::GetInstance()->m_dofs->GetIndex(node.first);
    vtkfile << stress(index, 1) << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;

  // output sigma_xy for two-dimension problem
  if (stress.cols() == 3)
  {
    vtkfile << "<DataArray type=\"Float64\" Name=\"sigma_xy\""
            << " NumberOfComponents=\"1\" format=\"ascii\" >" << std::endl;
    for (auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords)
    {
      int index = GlobalData::GetInstance()->m_dofs->GetIndex(node.first);
      vtkfile << stress(index, 2) << std::endl;
    }
    vtkfile << "</DataArray>" << std::endl;
  }
}

void MeshWriter::WritePointCoordsData(std::fstream &vtkfile)
{
  vtkfile << "<Points>" << std::endl;
  vtkfile << "<DataArray type=\"Float64\" Name=\"Points\""
          << " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  for(auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords)
  {
    for(int i = 0; i < node.second.size(); i++)
    {
      vtkfile << (node.second)(i) << "  ";
    }
    if(2 == node.second.size())
      vtkfile << "  0.";
    vtkfile << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;
  vtkfile << "</Points>" << std::endl;
}

void MeshWriter::WriteCellInfo(std::fstream &vtkfile)
{
  vtkfile << "<Cells>" << std::endl;
  
  // output element connectivity
  WriteElemConnecData(vtkfile);

  WriteElemOffsetData(vtkfile);

  WriteElemTypeData(vtkfile);

  vtkfile << "</Cells>" << std::endl;}

void MeshWriter::WriteElemConnecData(std::fstream &vtkfile)
{
  vtkfile << "<DataArray type=\"Int64\" Name=\"connectivity\""
          << " format=\"ascii\">" << std::endl;

  std::vector<int> elemNode;
  MatrixXd elemCoords;
  // std::string elemType = ShapeFunctions::GetElemType(elemDat->m_coords);
  for(auto elem : GlobalData::GetInstance()->m_elements->IterElementGroup(m_elementGroup))
  {
    elemNode = GlobalData::GetInstance()->m_elements->GetElementPtr()[elem]->GetNodes();
    elemCoords = GlobalData::GetInstance()->m_nodes->GetNodeCoords(elemNode);
    m_elemType.emplace_back( ShapeFunctions::GetElemType(elemCoords));
    
    for(auto nodeID : elemNode)
      vtkfile << GlobalData::GetInstance()->m_dofs->GetIndex(nodeID) << "  "; 

    vtkfile << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;
}

void MeshWriter::WriteElemOffsetData(std::fstream &vtkfile)
{
  vtkfile << "<DataArray type=\"Int64\" Name=\"offsets\""
          << " format=\"ascii\">" << std::endl;
  
  int count = 1;
  std::vector<int> elemNode;
  for(auto elem : GlobalData::GetInstance()->m_elements->IterElementGroup(m_elementGroup))
  {
    elemNode = GlobalData::GetInstance()->m_elements->GetElementPtr()[elem]->GetNodes();
    vtkfile << elemNode.size() * count << std::endl;
    count += 1;
  }
  vtkfile << "</DataArray>" << std::endl;
}

void MeshWriter::WriteElemTypeData(std::fstream &vtkfile)
{
  vtkfile << "<DataArray type=\"UInt8\" Name=\"types\""
          << " format=\"ascii\" RangeMin=\"9\" RangeMax=\"9\">" << std::endl;
  for(int i = 0; i < GlobalData::GetInstance()->m_elements->ElemGroupCount(m_elementGroup); i++)
    vtkfile << m_cellType.at(m_elemType[i]) << std::endl;
  vtkfile << "</DataArray>" << std::endl;
}

void MeshWriter::WritePVDFile()
{
  std::string filename = GlobalData::GetInstance()->m_prefix + ".pvd";
  std::fstream pvdfile;
  pvdfile.open(filename, std::ios::out);

  pvdfile << "<VTKFile byte_order=\'LittleEndian\' type=\'Collection\' version=\'0.1\'>" << std::endl;
  pvdfile << "<Collection>" << std::endl;
  
  for(int i = 0; i < m_k + 1; i++)
    pvdfile << "<DataSet file='" << GlobalData::GetInstance()->m_prefix
            << "-" << i << ".vtk ' groups='' part='0' timestep='" << i << "'/>" << std::endl;

  pvdfile << "</Collection>" << std::endl;
  pvdfile << "</VTKFile>" << std::endl;

  pvdfile.close();
}

void MeshWriter::WriteVTKFile()
{
  std::string filename = GlobalData::GetInstance()->m_prefix + "-" 
                       + std::to_string(m_k) + ".vtu";

  std::fstream vtkfile;
  vtkfile.open(filename, std::ios::out);
  
  vtkfile << "<?xml version=\"1.0\"?>" << std::endl;
  vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\""
          <<" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
  vtkfile << "<UnstructuredGrid>" << std::endl;
  
  vtkfile << "<Piece NumberOfPoints=\""
          << std::to_string(GlobalData::GetInstance()->m_nodes->m_nodeCoords.size())
          << "\" NumberOfCells=\""
          << std::to_string(GlobalData::GetInstance()->m_elements->ElemGroupCount(m_elementGroup))
          << "\">" << std::endl;

  WritePointInfo(vtkfile);

  WriteCellInfo(vtkfile);
  vtkfile << "</Piece>" << std::endl;
  vtkfile << "</UnstructuredGrid>" << std::endl;
  vtkfile << "</VTKFile>" << std::endl;
  
  vtkfile.close();
}