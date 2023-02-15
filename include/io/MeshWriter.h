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
  int m_k = 0;                          // Load Step
  int m_interval = 1;
  std::string m_fileName;               // Output File Name
  std::string m_elementGroup = "All";
};

ReflectRegister(MeshWriter, const nlohmann::json&)

#endif // MESHWRITER_H