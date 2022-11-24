#ifndef BASEMATERIAL_H
#define BASEMATERIAL_H

#include <vector>
#include <unordered_map>

#include <nlohmann/json.hpp>

class BaseMaterial
{
public:
  BaseMaterial(const nlohmann::json &props);
  virtual ~BaseMaterial();

protected:
  void SetIter(int iIter);

  int SetHistoryParameter(const std::string &name, double vale);

  double GetHistoryParameter(const std::string&name);

  void CommitHistory();

  virtual void ComputeDMatrix() = 0;

protected:
  int m_iIter = -1;
  std::vector<std::unordered_map<std::string, double>> m_current;
  std::vector<std::unordered_map<std::string, double>> m_history;
  std::unordered_map<std::string, double> m_initHistory;
  nlohmann::json m_props;
};

#endif // BASEMATERIAL_H