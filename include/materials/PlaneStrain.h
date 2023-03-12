#ifndef PLEANESTRAIN_H
#define PLEANESTRAIN_H

#include <materials/BaseMaterial.h>
#include <util/ObjectFactory.h> 

class PlaneStrain : public BaseMaterial
{
public:
  /**
   * @Brief: Construct a new Plane Strain object
   * 
   * @param props 
   */
  PlaneStrain(const nlohmann::json &props);

  /**
   * @Brief: Destroy the Plane Strain object
   * 
   */
  ~PlaneStrain();

private:
  virtual void ComputeDMatrix() override;

private:
  double m_E = 0.;
  double m_nu = 0.;
};
ReflectRegister(PlaneStrain, const nlohmann::json &)

#endif // PLEANESTRAIN_H