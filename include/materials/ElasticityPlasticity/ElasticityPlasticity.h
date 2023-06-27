#ifndef ELASTICITYPLASTICITY_H
#define ELASTICITYPLASTICITY_H

#include <materials/BaseMaterial.h>
#include <util/ObjectFactory.h>
#include <materials/LinearElasticity.h>
#include <util/Math.h>

class ElasticityPlasticity : public BaseMaterial
{
public:
  ElasticityPlasticity(const nlohmann::json &matProps);

  ~ElasticityPlasticity();

  virtual VectorXd GetStress(const std::shared_ptr<Kinematics> &kin,
                             const VectorXd &increDisp,
                             const MatrixXd &dphi = MatrixXd::Zero(0, 0)) override;

private:
  virtual void ComputeDMatrix() override;

  void StressRotation(VectorXd &stress,
                      const VectorXd &increDisp,
                      const MatrixXd &dphi);

private:
  double m_K;
  double m_G;
  double m_plaMod;
  double m_yieldStress;
  double m_waveSpeed;
  std::string m_rateType;
  VectorXd m_oneVec;

  std::shared_ptr<LinearElasticity> lineMat;
};
ReflectRegister(ElasticityPlasticity, const nlohmann::json &)

#endif // ELASTICITYPLASTICITY_H