#ifndef GRAPHWRITER_H
#define GRAPHWRITER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>
#include <petscvec.h>
#include <fstream>

class GraphWriter : public BaseModule
{
public:
  GraphWriter(const nlohmann::json &props);
  ~GraphWriter();

  virtual void Run() override;

private:
  Vec& GetGlobalData(const std::string &name);

private:
  std::vector<std::string> m_columns;
  std::vector<std::vector<double>> m_data;
   std::fstream m_outFileStream;
  
  int index;
  double tempDouble;
  double factor;
  std::string dofType;
};
ReflectRegister(GraphWriter, const nlohmann::json &);

#endif // GRAPHWRITER_H