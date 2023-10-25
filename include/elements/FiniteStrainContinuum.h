/**
 * @File Name:     FiniteStrainContinuum.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef FINITESTRAINCONTINUUM_H
#define FINITESTRAINCONTINUUM_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

/**
 * @Brief:  Total Lagrange Formulation
 * 
 */

class FiniteStrainContinuum : public Element
{
public:
  /**
   * @Brief: Construct a new Finite Strain Continuum object
   * 
   * @param elemShape 
   * @param elemNodes 
   * @param modelProps 
   */
  FiniteStrainContinuum(const std::string &elemShape,
                        const std::vector<int> &elemNodes,
                        const nlohmann::json &modelProps);
  
  /**
   * @Brief:         Destroy the Finite Strain Continuum object
   * 
   */
  ~FiniteStrainContinuum();

  /**
   * @Brief: Compute the Element Tangent Stiffness
   * 
   * @param elemDat 
   */
  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  /**
   * @Brief: Compute the Kinematics Variable Such as Deformation Gradient,
   *          Cauchy-Green Strain Tensor and Strain Vector
   * 
   * @param dphi 
   * @param elState 
   */
  void GetKinematics(const MatrixXd&dphi,
                     const VectorXd&elState);
  
  /**
   * @Brief: Compute Strain Matrix B
   * 
   * @param dphi 
   * @param F 
   */
  void GetBMatrix(const MatrixXd&dphi,
                  const MatrixXd&F);
  
  /**
   * @Brief: Form Stress Matrix to Compute Nonlinear Stiffness Matrix
   * 
   * @param stress 
   */
  void Stress2Matrix(const VectorXd&stress);

  /**
   * @Brief: Compute Nonlinear Strain Matrix BNL
   * 
   * @param dphi 
   */
  void GetBNLMatrix(const MatrixXd&dphi);
  
  /**
   * @Brief: Initialize Some Basic Variables
   * 
   * @param elemShape 
   */
  void Initialize(const std::string &elemShape);

  // virtual void SetElemNodeOrdered();

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  
  double detJac = 0.;                                   // the Determinant of Jacobian Matrix
  /*jac = [pXpxi1 pXpxi2 pXpxi3,
           pYpxi1 pYpxi2 pYpxi3,
           pZpxi1 pZpxi2 pZpxi3];*/
  MatrixXd jac;                                         // the Jacobian Matrix
  MatrixXd invJac;                                      // the Inverse Jacobian Matrix

  /* pHpX = [pH1pX pH1pY pH1pZ,
             pH2pX pH2pY pH2pZ,
             ...
             pHnpX pHnpY pHnpZ]*/
  MatrixXd pHpX;                                        // the Derivative of Shape Function
  MatrixXd B;                                           // the Strain Matrix
  VectorXd sigma;                                       // the Stress Vector
  MatrixXd D;                                           // the Tangent Matrix
  MatrixXd Bnl;

  // stress matrix
  MatrixXd T; 
};

ReflectRegister(FiniteStrainContinuum, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUM_H