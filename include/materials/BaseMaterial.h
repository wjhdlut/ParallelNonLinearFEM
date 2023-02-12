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

  virtual std::vector<double> GetStress(const std::shared_ptr<Kinematics>&kin) = 0;

  void SetIter(int iIter);

  void CommitHistory();

  inline std::vector<std::vector<double>> GetTangMatrix(){
    return m_D;
  }

protected:
  int SetHistoryParameter(const std::string &name, double vale);

  double GetHistoryParameter(const std::string&name);

  virtual void ComputeDMatrix() = 0;

protected:
  int m_iIter = -1;
  std::vector<std::unordered_map<std::string, double>> m_current;
  std::vector<std::unordered_map<std::string, double>> m_history;
  std::unordered_map<std::string, double> m_initHistory;
  nlohmann::json m_props;
  std::vector<std::vector<double>> m_D;
};

#endif // BASEMATERIAL_H