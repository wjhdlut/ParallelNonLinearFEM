#ifndef BASEMATERIAL_H
#define BASEMATERIAL_H

#include <vector>
#include <unordered_map>

#include <nlohmann/json.hpp>
#include <util/Kinematics.h>

class BaseMaterial
{
public:
  BaseMaterial(const nlohmann::json &props);
  virtual ~BaseMaterial();

  virtual std::vector<double> GetStress(const std::shared_ptr<Kinematics>&kin,
                                        const std::vector<double> &increDisp = std::vector<double>(),
                                        const Matrix &dhpi = Matrix());

  void SetIter(int iIter);

  void CommitHistory();

  double SetMaterialParamter(const std::string &name);

  inline std::vector<std::vector<double>> GetTangMatrix(){
    if(0 == m_D.size()) ComputeDMatrix();
    return m_D;
  }

  inline double ReturnRho()
  {
    return m_rho;
  }

protected:
  int SetHistoryParameter(const std::string &name, double value);

  double GetHistoryParameter(const std::string&name);

  virtual void ComputeDMatrix() = 0;

protected:
  int m_iIter = -1;
  double m_rho = 0.;
  double m_E = 0.;
  double m_nu = 0.;
  std::vector<std::unordered_map<std::string, double>> m_current;
  std::vector<std::unordered_map<std::string, double>> m_history;
  std::unordered_map<std::string, double> m_initHistory;
  nlohmann::json m_props;
  std::vector<std::vector<double>> m_D;
};

#endif // BASEMATERIAL_H