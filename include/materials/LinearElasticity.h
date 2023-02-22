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
  void ComputeDMatrix();

private:
  double m_E = 0.;
  double m_nu = 0.;
};

ReflectRegister(LinearElasticity, const nlohmann::json &)

#endif // LINEARELASTICITY_H