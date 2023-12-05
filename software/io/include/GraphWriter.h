/**
 * @File Name:     GraphWriter.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef GRAPHWRITER_H
#define GRAPHWRITER_H

#include <petscvec.h>
#include <fstream>
#include <eigen3/Eigen/Dense>

#include "../../util/include/BaseModule.h"
#include "../../util/include/ObjectFactory.h"

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