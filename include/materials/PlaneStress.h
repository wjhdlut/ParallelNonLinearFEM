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

  /**
   * @Brief: Compute the Stress Vector
   * 
   * @param kin 
   * @return std::vector<double> 
   */
  virtual std::vector<double> GetStress(const std::shared_ptr<Kinematics>&kin) override;

private:
  /**
   * @Brief:  Compute the Tangent Modulue Matrix
   * 
   */
  void ComputeDMatrix() override;

private:
  double m_E = 0.;
  double m_nu = 0.;
};

ReflectRegister(PlaneStress, const nlohmann::json &)
#endif // PLANESTRESS_H