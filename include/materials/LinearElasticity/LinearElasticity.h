#ifndef LINEARELASTICITY_H
#define LINEARELASTICITY_H


#include <materials/BaseMaterial.h>
#include <util/ObjectFactory.h>

class LinearElasticity : public BaseMaterial
{
public:
  /**
   * @Brief: Construct a new Linear Elasticity object
   * 
   * @param props 
   */
  LinearElasticity(const nlohmann::json &props);

  /**
   * @Brief: Destroy the Linear Elasticity object
   * 
   */
  ~LinearElasticity();

private:
  /**
   * @Brief:  Compute the Tangent Modulue Matrix
   * 
   */
  virtual void ComputeDMatrix() override;

private:
  /**
   * @Brief: For 3D Problem
   * 
   */
  void For3D();
  
  /**
   * @Brief: For Plane Strain Problem
   * 
   */
  void ForPlaneStrain();

  /**
   * @Brief: For Plane Stress Problem
   * 
   */
  void ForPlaneStress();
};
ReflectRegister(LinearElasticity, const nlohmann::json &)

#endif // LINEARELASTICITY_H