
#ifndef SMALLSTRAINCONTINUUM_H
#define SMALLSTRAINCONTINUUM_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class SmallStrainContinuum : public Element
{
public:
  SmallStrainContinuum(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  
  ~SmallStrainContinuum();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

  virtual void GetMassMatrix(std::shared_ptr<ElementData> &elemDat) override {}

private:
  std::shared_ptr<Kinematics> GetKinematics(const Matrix &B, const std::vector<double> &elState);

  Matrix GetBMatrix(const Matrix &dphi);
};

ReflectRegister(SmallStrainContinuum, const std::vector<int> &, const nlohmann::json &)

#endif // SMALLSTRAINCONTINUUM_H