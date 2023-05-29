#ifndef MESHWRITER_H
#define MESHWRITER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>

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
  void WriteVTKFile();

  void WritePVDFile();

  void WritePointInfo(std::fstream &vtkfile);

  void WriteCellInfo(std::fstream &vtkfile);

  void WritePointDisplaceData(std::fstream &vtkfile);

  void WritePointStressData(std::fstream &vtkfile);

  void WritePointCoordsData(std::fstream &vtkfile);

  void WriteElemConnecData(std::fstream &vtkfile);

  void WriteElemOffsetData(std::fstream &vtkfile);

  void WriteElemTypeData(std::fstream &vtkfile);

private:
  int m_k = 0;                          // Load Step
  int m_interval = 1;
  std::string m_fileName;               // Output File Name
  std::string m_elementGroup = "All";
  std::vector<std::string> m_elemType;

  const std::unordered_map<std::string, std::string> m_cellType = {{"Hexa8", "12"},
                                                                   {"Quad4", "9"},
                                                                   {"Quad8", "9"}};
};

ReflectRegister(MeshWriter, const nlohmann::json&)

#endif // MESHWRITER_H