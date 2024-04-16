/**
 * @File Name:     MeshWriter.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef MESHWRITER_H
#define MESHWRITER_H

#include "../../util/include/BaseModule.h"
#include "../../util/include/ObjectFactory.h"

class MeshWriter : public BaseModule
{
public:
  /**
   * @Brief: Construct Output Mesh object
   * 
   * @param props 
   */
  MeshWriter(const nlohmann::json&props);

  /**
   * @Brief: Destroy Output Mesh object
   * 
   */
  ~MeshWriter();

  /**
   * @Brief:  Output Mesh Information
   * 
   */
  virtual void Run() override;

private:
  /**
   * @Brief: Write Resulte into VTK Files
   * 
   */
  void WriteVTKFile();

  /**
   * @Brief: Write a PVD File
   * 
   */
  void WritePVDFile();

  /**
   * @Brief: Write Nodal Displacement and Stress
   * 
   * @param vtkfile 
   */
  void WritePointInfo(std::fstream &vtkfile);

  /**
   * @Brief: Write Elemental Information
   * 
   * @param vtkfile 
   */
  void WriteCellInfo(std::fstream &vtkfile);

  /**
   * @Brief: Write Nodal Displacement Result
   * 
   * @param vtkfile 
   */
  void WritePointDisplaceData(std::fstream &vtkfile);

  /**
   * @Brief: Write Nodal Stress Reault
   * 
   * @param vtkfile 
   */
  void WritePointStressData(std::fstream &vtkfile);

  /**
   * @Brief: Write Nodal Coordinates Reslut
   * 
   * @param vtkfile 
   */
  void WritePointCoordsData(std::fstream &vtkfile);

  /**
   * @Brief: Write Element Connection Information
   * 
   * @param vtkfile 
   */
  void WriteElemConnecData(std::fstream &vtkfile);

  /**
   * @Brief: Write Element Offset Data
   * 
   * @param vtkfile 
   */
  void WriteElemOffsetData(std::fstream &vtkfile);

  /**
   * @Brief: Write Element Type Information
   * 
   * @param vtkfile 
   */
  void WriteElemTypeData(std::fstream &vtkfile);

private:
  int m_k = 0;                                           // Load Step
  int m_interval = 1;
  std::string m_fileName;                                // Output File Name
  std::string m_elementGroup = "All";
  std::vector<std::string> m_elemType;

  const std::unordered_map<std::string, std::string> m_cellType = {{"Hexa8", "12"},
                                                                   {"Quad4", "9"},
                                                                   {"Quad8", "9"},
                                                                   {"Tria3", "9"}};
};

ReflectRegister(MeshWriter, const nlohmann::json&)

#endif // MESHWRITER_H