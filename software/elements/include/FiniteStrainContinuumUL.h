/**
 * @File Name:     FiniteStrainContinuumUL.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef FINITESTRAINCONTINUUMUL_H
#define FINITESTRAINCONTINUUMUL_H

#include "Element.h"
#include "../../util/include/ObjectFactory.h"

/**
 * @Brief: Update Lagrange Formulation
 * 
 */

class FiniteStrainContinuumUL : public Element
{
public:
  /**
   * @Brief: Construct a new Finite Strain Continuum Element
   *         Base on UL Formulation
   * 
   * @param elemShape 
   * @param elemNodes 
   * @param modelProps 
   */
  FiniteStrainContinuumUL(const std::string &elemShape,
                          const std::vector<int> &elemNodes,
                          const nlohmann::json &modelProps);

  /**
   * @Brief: Destroy the Finite Strain Continuum Element
   *         Based on UL Formulation
   * 
   */
  ~FiniteStrainContinuumUL();

  /**
   * @Brief: Compute the Tangent Stiffness
   * 
   * @param elemDat 
   */
  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  /**
   * @Brief: Compute the Kinematics Variables, Such as Deformation Gradient,
   *         Cauchy-Green Strain Tensor and Strain Vector
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
   */
  void GetBMatrix(const MatrixXd&dphi);
  
  /**
   * @Brief: Form the Stress Matrix to Compute Nonlinear Stiffness Matrix
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
   * @Brief:  Update the Element Node Coordinate
   * 
   * @param elemDat 
   * @return MatrixXd 
   */
  MatrixXd UpdateNodeCoords(std::shared_ptr<ElementData>&elemDat);

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  
  /*jac = [pxpxi1 pxpxi2 pxpxi3,
           pypxi1 pypxi2 pypxi3,
           pzpxi1 pzpxi2 pzpxi3];*/
  MatrixXd jac;                                         // the Jacobian Matrix

  /*jac = [pXpxi1 pXpxi2 pXpxi3,
           pYpxi1 pYpxi2 pYpxi3,
           pZpxi1 pZpxi2 pZpxi3];*/
  MatrixXd Jac;                                         // the Jacobian Matrix
  
  MatrixXd invJac;                                      // the Inverse Jacobian Matrix

  /* pHpX = [pH1pX pH1pY pH1pZ,
             pH2pX pH2pY pH2pZ,
             ...
             pHnpX pHnpY pHnpZ]*/
  MatrixXd pHpX;

  /* pHpX = [pH1px pH1py pH1pz,
             pH2px pH2py pH2pz,
             ...
             pHnpx pHnpy pHnpz]*/
  MatrixXd pHpx;                                        // the Derivative of Shape Function
  MatrixXd B;                                           // the Strain Matrix
  VectorXd sigma;                                       // the Stress Vector
  MatrixXd D;                                           // the Tangent Matrix
  MatrixXd Bnl;
  MatrixXd T;                                           // stress matrix
  MatrixXd nodeCoord;
};

ReflectRegister(FiniteStrainContinuumUL, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUMUL_H