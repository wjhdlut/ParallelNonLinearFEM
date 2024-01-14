/**
 * @File Name:     SmallStrainContinuum.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef SMALLSTRAINCONTINUUM_H
#define SMALLSTRAINCONTINUUM_H

#include "Element.h"
#include "../../util/include/ObjectFactory.h"

#include <fstream>

class SmallStrainContinuum : public Element
{
public:
  /**
   * @Brief: Construct a new Small Strain Continuum object
   * 
   * @param elemShape 
   * @param elemNodes 
   * @param modelProps 
   */
  SmallStrainContinuum(const std::string &elemShape,
                       const std::vector<int> &elemNodes,
                       const nlohmann::json &modelProps);
  
  /**
   * @Brief: Destroy the Small Strain Continuum object
   * 
   */
  ~SmallStrainContinuum();

  /**
   * @Brief: Compute the Tangent Stiffness Matrix
   * 
   * @param elemDat 
   */
  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:  
  /**
   * @Brief: Compute the Kinematics Variables, Such as Deformation Gradient,
   *         Cauchy-Green Strain Tensor and Strain Vector
   * 
   * @param elState 
   */
  void GetKinematics(const std::shared_ptr<ElementData> &elemDat);

  /**
   * @Brief: Compute Strain Matrix B
   * 
   * @param dphi 
   */
  void GetBMatrix(const MatrixXd &dphi);

  /**
   * @Brief: Initialize Some Variables
   * 
   * @param elemShape 
   */
  void Initialize(const std::string &elemShape);

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  
  double detJac = 0.;                                   // the Determinant of Jacobian Matrix
  /*jac = [pXpxi1 pXpxi2 pXpxi3,
           pYpxi1 pYpxi2 pYpxi3,
           pZpxi1 pZpxi2 pZpxi3];*/
  MatrixXd jac;                                        // the Jacobian Matrix
  MatrixXd invJac;                                     // the Inverse Jacobian Matrix

  /* pHpX = [pH1pX1 pH1pX2 pH1pX3,
             pH2pX1 pH2pX2 pH2pX3,
             ...
             pHnpX1 pHnpX2 pHnpX3,]*/
  MatrixXd pHpX;                                       // the Derivative of Shape Function
  MatrixXd B;                                          // the Strain Matrix
  VectorXd sigma;                                      // the Stress Vector
  MatrixXd D;                                          // the Tangent Matrix
};

ReflectRegister(SmallStrainContinuum, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif // SMALLSTRAINCONTINUUM_H