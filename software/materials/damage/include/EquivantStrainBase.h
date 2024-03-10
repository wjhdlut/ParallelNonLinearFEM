/**
 * @File Name:     EquivanceBase.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2024-01-23
 * 
 * @Copyright Copyright (c) 2024 JianHuaWang
 * 
 */

#ifndef EQUIVANTSTRAINBASE_H
#define EQUIVANTSTRAINBASE_H

#include "../../../nlohmann/json.hpp"
#include "../../../util/include/Kinematics.h"

class EquivantStrainBase
{
public:
  /**
   * @Brief: Default Constructor of EquivantStrainBase Class
   * 
   */
  EquivantStrainBase() = default;

  /**
   * @Brief: Constructor of EquivantStrainBase Class
   * 
   * @param props 
   */
  EquivantStrainBase(const nlohmann::json&props);

  /**
   * @Brief: Destructor of EquivantStrainBase Class
   * 
   */
  ~EquivantStrainBase();

public:
  /**
   * @Brief: Compute Equivance Strain
   * 
   * @param kin 
   */
  virtual void CompEquivanceStrain(double &eps,
                                   VectorXd &dEpsdStrain,
                                   const std::shared_ptr<Kinematics>&kin) = 0;

  /**
   * @Brief: Read Some Basic Variables
   * 
   */
  virtual void ReadParameters() = 0;

private:
  /**
   * @Brief: Initialize Some Basic Variables
   * 
   */
  void Initialize();

protected:
  bool m_planeStrainFlag = false;
  bool m_planeStressFlag = false;
  nlohmann::json m_props;
};


#endif //EQUIVANTSTRAINBASE_H