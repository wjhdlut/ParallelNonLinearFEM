/**
 * @File Name:     VonMisesEquivantStrain.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2024-01-23
 * 
 * @Copyright Copyright (c) 2024 JianHuaWang
 * 
 */

/*********************************************************************
 * Compute Equivance Strain Based on Modified Von Mises Definition
 * Referencing Book [E.q. (6.21) in Chapter 6.2]：
 *  'Non-Linear Finite Element Analysis of Solids and Structures'
 *   R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
 *   John Wiley and Sons, 2012, ISBN 978-0470666449  
 *********************************************************************/

#ifndef VONMISESEQUIVANTSTRAIN_H
#define VONMISESEQUIVANTSTRAIN_H

#include "EquivantStrainBase.h"
#include "../../../util/include/ObjectFactory.h"

class VonMisesEquivantStrain : public EquivantStrainBase
{
public:
  /**
   * @Brief: Default Constructor of VonMisesEquivantStrain Class
   * 
   */
  VonMisesEquivantStrain() = default;

  /**
   * @Brief: Constructor of VonMisesEquivantStrain Class
   * 
   * @param props 
   */
  VonMisesEquivantStrain(const nlohmann::json &props);

  /**
   * @Brief: Destructor of VonMisesEquivantStrain Class
   * 
   */
  ~VonMisesEquivantStrain();

public:
  /**
   * @Brief: Compute the Equivance Strain
   * 
   * @param kin 
   */
  virtual void CompEquivanceStrain(double &eps,
                                   VectorXd &dEpsdStrain,
                                   const std::shared_ptr<Kinematics>&kin);

  /**
   * @Brief: Read SOme Basic Variables
   * 
   */
  virtual void ReadParameters();

private:
  /**
   * @Brief: Equivance Strain for Plane Strain Problem
   * 
   * @param kin 
   */
  Vector3d PlaneStrainEquivanceStrain(double &eps,
                                      const std::shared_ptr<Kinematics>&kin);

  /**
   * @Brief: Equivance Strain for Plane Stress Problem
   * 
   * @param kin 
   */
  Vector3d PlaneStressEquivandeStrain(double &eps,
                                      const std::shared_ptr<Kinematics>&kin);

  /**
   * @Brief: Equivance Strain for Three Dimensional Problem
   * 
   * @param kin 
   */
  VectorXd ThreeDimensionEquivandeStrain(double &eps,
                                         const std::shared_ptr<Kinematics>&kin);

  /**
   * @Brief: Initialize Some Basic Variables
   * 
   */
  void Initialize();

private: 
  double m_nu = 0.; 
  double m_a1 = 0.;
  double m_a2 = 0.;
  double m_a3 = 0.;
  double m_a4 = 0.;
  double m_c  = 0;
  double m_sc = 1./3.;
  double m_k  = 0.;
};

ReflectRegister(VonMisesEquivantStrain, const nlohmann::json&)

#endif // VONMISESEQUIVANTSTRAIN_H