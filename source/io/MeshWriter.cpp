
#include <fstream>

#include <io/MeshWriter.h>
#include <util/DataStructure.h>
#include <iostream>

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

  vtkfile << "<PointData>" << std::endl;
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
    if(node.second.size() < 3) vtkfile << "  0." << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;

  // output stress information
  // output sigma_xx
  Matrix stress = GlobalData::GetInstance()->GetData("stresses");
  vtkfile << "<DataArray type=\"Float64\" Name=\"sigma_xx\""
          << " NumberOfComponents=\"1\" format=\"ascii\" >" << std::endl;
  for (auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords)
  {
    int index = GlobalData::GetInstance()->m_dofs->GetIndex(node.first);
    vtkfile << stress[index][0] << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;

  // output sigma_yy
  vtkfile << "<DataArray type=\"Float64\" Name=\"sigma_yy\""
          << " NumberOfComponents=\"1\" format=\"ascii\" >" << std::endl;
  for (auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords)
  {
    int index = GlobalData::GetInstance()->m_dofs->GetIndex(node.first);
    vtkfile << stress[index][1] << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;

  // output sigma_xy for two-dimension problem
  if (stress[0].size() == 3)
  {
    vtkfile << "<DataArray type=\"Float64\" Name=\"sigma_xy\""
            << " NumberOfComponents=\"1\" format=\"ascii\" >" << std::endl;
    for (auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords)
    {
      int index = GlobalData::GetInstance()->m_dofs->GetIndex(node.first);
      vtkfile << stress[index][2] << std::endl;
    }
    vtkfile << "</DataArray>" << std::endl;
  }

  vtkfile << "</PointData>" << std::endl;
  vtkfile << "<CellData>" << std::endl;
  vtkfile << "</CellData>" << std::endl;

  // output node coordinate information
  vtkfile << "<Points>" << std::endl;
  vtkfile << "<DataArray type=\"Float64\" Name=\"Points\""
          << " NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  for(auto node : GlobalData::GetInstance()->m_nodes->m_nodeCoords)
  {
    for(auto coord : node.second)
    {
      vtkfile << coord << "  ";
    }
    if(2 == node.second.size())
      vtkfile << "  0.";
    vtkfile << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;
  vtkfile << "</Points>" << std::endl;
  vtkfile << "<Cells>" << std::endl;
  
  // output element connectivity
  vtkfile << "<DataArray type=\"Int64\" Name=\"connectivity\""
          << " format=\"ascii\">" << std::endl;

  std::vector<int> elemNode;
  for(auto elem : GlobalData::GetInstance()->m_elements->IterElementGroup(m_elementGroup))
  {
    elemNode = GlobalData::GetInstance()->m_elements->GetElementPtr()[elem]->GetNodes();

    for(auto nodeID : elemNode)
      vtkfile << GlobalData::GetInstance()->m_dofs->GetIndex(nodeID) << "  "; 

    vtkfile << std::endl;
  }
  vtkfile << "</DataArray>" << std::endl;

  vtkfile << "<DataArray type=\"Int64\" Name=\"offsets\""
          << " format=\"ascii\">" << std::endl;
  int count = 1;
  for(auto elem : GlobalData::GetInstance()->m_elements->IterElementGroup(m_elementGroup))
  {
    elemNode = GlobalData::GetInstance()->m_elements->GetElementPtr()[elem]->GetNodes();
    vtkfile << elemNode.size() * count << std::endl;
    count += 1;
  }
  vtkfile << "</DataArray>" << std::endl;

  vtkfile << "<DataArray type=\"UInt8\" Name=\"types\""
          << " format=\"ascii\" RangeMin=\"9\" RangeMax=\"9\">" << std::endl;
  for(int i = 0; i < GlobalData::GetInstance()->m_elements->ElemGroupCount(m_elementGroup); i++)
    vtkfile << 9 << std::endl;
  vtkfile << "</DataArray>" << std::endl;

  vtkfile << "</Cells>" << std::endl;
  vtkfile << "</Piece>" << std::endl;
  vtkfile << "</UnstructuredGrid>" << std::endl;
  vtkfile << "</VTKFile>" << std::endl;
  
  vtkfile.close();

  // write pvd file
  filename = GlobalData::GetInstance()->m_prefix + ".pvd";
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

  m_k += 1;
}