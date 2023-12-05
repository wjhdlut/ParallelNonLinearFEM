#ifndef PLANESTRESS_H
#define PLANESTRESS_H

#include <materials/BaseMaterial.h>
#include <util/ObjectFactory.h>

class PlaneStress : public BaseMaterial
{
public:
  /**
   * @Brief: Construct a new Plane Stress Materials Object
   * 
   * @param props 
   */
  PlaneStress(const nlohmann::json &props);

  /**
   * @Brief: Destroy the Plane Stress Material Object
   * 
   */
  ~PlaneStress();

private:
  /**
   * @Brief:  Compute the Tangent Modulue Matrix
   * 
   */
  virtual void ComputeDMatrix() override;
};

ReflectRegister(PlaneStress, const nlohmann::json &)
#endif // PLANESTRESS_H