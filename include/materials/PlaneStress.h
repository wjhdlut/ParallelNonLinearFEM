#ifndef PLANESTRESS_H
#define PLANESTRESS_H

#include <materials/BaseMaterial.h>
#include <util/ObjectFactory.h>

class PlaneStress : public BaseMaterial
{
public:
  PlaneStress(const nlohmann::json &props);
  ~PlaneStress();

private:
  void ComputeDMatrix() override;

private:
  double m_E = 0.;
  double m_nu = 0.;
  double m_D[3][3] = {{0.}};
};

ReflectRegister(PlaneStress, const nlohmann::json &)
#endif // PLANESTRESS_H