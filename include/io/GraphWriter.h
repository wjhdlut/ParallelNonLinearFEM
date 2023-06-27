#ifndef GRAPHWRITER_H
#define GRAPHWRITER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>
#include <petscvec.h>
#include <fstream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

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
  std::vector<VectorXd> m_data;
   std::fstream m_outFileStream;
  
  int index;
  double tempDouble;
  double factor;
  std::string dofType;
};
ReflectRegister(GraphWriter, const nlohmann::json &);

#endif // GRAPHWRITER_H