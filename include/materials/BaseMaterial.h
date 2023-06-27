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

  virtual VectorXd GetStress(const std::shared_ptr<Kinematics>&kin,
                             const VectorXd &increDisp = VectorXd::Zero(0),
                             const MatrixXd &dhpi = MatrixXd::Zero(0, 0));

  void SetIter(int iIter);

  void CommitHistory();

  double SetMaterialParamter(const std::string &name);

  inline MatrixXd GetTangMatrix(){
    if(0 == m_D.rows()) ComputeDMatrix();
    return m_D;
  }
  
  inline MatrixXd GetTangentMatrix(const MatrixXd &F)
  {
    if(0 == m_D.rows()) ComputeDMatrix();
    TangentDMatrixToInitial(F);
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

  void TangentDMatrixToCurrent(const MatrixXd &F);

  void TangentDMatrixToInitial(const MatrixXd &F);

private:
  void GetTransMatrix(MatrixXd &T, const MatrixXd &F);

protected:
  int m_iIter  = -1;
  double m_rho = 0.;
  double m_E   = 0.;
  double m_nu  = 0.;
  std::vector<std::unordered_map<std::string, double>> m_current;
  std::vector<std::unordered_map<std::string, double>> m_history;
  std::unordered_map<std::string, double> m_initHistory;
  nlohmann::json m_props;
  MatrixXd m_D;
};

#endif // BASEMATERIAL_H