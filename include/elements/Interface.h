/**
 * @File Name:     Interface.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef INTERFACE_H
#define INTERFACE_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class Interface : public Element
{
public:
  /**
   * @Brief: Construct a new Interface object
   * 
   * @param elemShape 
   * @param elemNodes 
   * @param modelProps 
   */
  Interface(const std::string &elemShape,
            const std::vector<int> &elemNodes,
            const nlohmann::json &modelProps);
  
  /**
   * @Brief: Destroy the Interface object
   * 
   */
  ~Interface();

  /**
   * @Brief: Compute the Tangent Stiffness Matrix
   * 
   * @param elemDat 
   */
  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  /**
   * @Brief: Compute the Rotation Matrix
   * 
   * @param coords 
   * @param state 
   * @return MatrixXd 
   */
  MatrixXd GetRotation(const MatrixXd &coords, const VectorXd &state);
  
  /**
   * @Brief: Compute Strain Matrix B
   * 
   * @param H 
   * @param R 
   */
  void GetBMatrix(const VectorXd &H, const MatrixXd &R);

  /**
   * @Brief: Compute the Kinematics Variable Such as Deformation Gradient,
   *         Cauchy-Green Strain Tensor and Strain Vector
   * 
   * @param elState 
   */
  void GetKinematics(const VectorXd &elState);

  /**
   * @Brief: Initialize Some Variables of Interface Element
   * 
   * @param elemShape 
   */
  void Initialize(const std::string &elemShape);

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  MatrixXd B;                                 // the Strain Matrix
  VectorXd sigma;                             // the Stress Vector
  MatrixXd D;                                 // the Tangent Matrix
};

ReflectRegister(Interface, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif // INTERFACE_H