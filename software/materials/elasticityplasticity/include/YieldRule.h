/**
 * @File Name:     YieldRule.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-11-13
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef YIELDRULE_H
#define YIELDRULE_H

#include <eigen3/Eigen/Dense>

#include "../../../nlohmann/json.hpp"

using namespace Eigen;

class YieldRule
{
public:
  /**
   * @Brief: Construct a new Yield Rule object
   * 
   */
  YieldRule() {};
  
  /**
   * @Brief: Construct a new Yield Rule object
   * 
   * @param matProps 
   */
  YieldRule(const nlohmann::json &matProps){};

  /**
   * @Brief: Destroy the Yield Rule object
   * 
   */
  virtual ~YieldRule() {};
  
  /**
   * @Brief: Compute Yield Function value
   * 
   */
  virtual double CompYieldFunction(const VectorXd &stress) = 0;

protected:
  bool m_plainStrain = false;
  bool m_plainStress = false;
  double m_nu = 0.;
};

#endif