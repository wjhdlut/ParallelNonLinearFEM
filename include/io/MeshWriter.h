#ifndef MESHWRITER_H
#define MESHWRITER_H

#include <util/BaseModule.h>

class MeshWriter : public BaseModule
{
public:
  MeshWriter(const nlohmann::json&props);
  ~MeshWriter();

  void Run();

private:
  int m_k = 0;
  int m_interval = 0;
  std::string m_fileName;
  std::string m_elementGroup = "All";
};

#endif // MESHWRITER_H