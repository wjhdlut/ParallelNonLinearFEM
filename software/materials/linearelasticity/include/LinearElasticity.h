/**
 * @File Name:     LinearElasticity.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef LINEARELASTICITY_H
#define LINEARELASTICITY_H


#include "../../include/BaseMaterial.h"
#include "../../../util/include/ObjectFactory.h"

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

  /**
   * @Brief: For AxiSYmmetric Problem
   * 
   */
  void ForAxiSymmetry();
};
ReflectRegister(LinearElasticity, const nlohmann::json &)

#endif // LINEARELASTICITY_H