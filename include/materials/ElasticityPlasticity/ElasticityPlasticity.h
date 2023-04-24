#ifndef ELASTICITYPLASTICITY_H
#define ELASTICITYPLASTICITY_H

#include <materials/BaseMaterial.h>
#include <util/ObjectFactory.h>
#include <materials/LinearElasticity.h>

class ElasticityPlasticity : public BaseMaterial
{
public:
  ElasticityPlasticity(const nlohmann::json &matProps);

  ~ElasticityPlasticity();

  virtual std::vector<double> GetStress(const std::shared_ptr<Kinematics> &kin,
                                        const std::vector<double> &increDisp,
                                        const Matrix &dphi = Matrix()) override;

private:
  virtual void ComputeDMatrix() override;

  void StressRotation(std::vector<double> &stress,
                      const std::vector<double> &increDisp,
                      const Matrix &dphi);

private:
  double m_K;
  double m_G;
  double m_plaMod;
  double m_yieldStress;
  double m_waveSpeed;
  std::string m_rateType;
  std::vector<double> m_oneVec = {1., 1., 1., 0., 0., 0.};

  std::shared_ptr<LinearElasticity> lineMat;
};
ReflectRegister(ElasticityPlasticity, const nlohmann::json &)

#endif // ELASTICITYPLASTICITY_H