/**
 * @File Name:     KirchhoffBeam.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-03-13
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef KIRCHHOFFBEAM_H
#define KIRCHHOFFBEAM_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class KirchhoffBeam : public Element
{
public:
  /**
   * @Brief: Construct a new Kirchhoff Beam object
   * 
   * @param elemShape 
   * @param elemNode 
   * @param modelProps 
   */
  KirchhoffBeam(const std::string &elemShape,
                const std::vector<int> &elemNode,
                const nlohmann::json &modelProps);
  
  /**
   * @Brief: Destroy the Kirchhoff Beam object
   * 
   */
  ~KirchhoffBeam();

  /**
   * @Brief: Compute the Tangent Stiffness
   * 
   * @param elemDat 
   */
  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  /**
   * @Brief: Transform Coordinate into Element Coordinate System
   * 
   * @param a 
   * @param coord 
   * @return VectorXd 
   */
  VectorXd TolemCoordinate(const VectorXd &a, const MatrixXd &coord);

  /**
   * @Brief: Compute the Bu Vector
   * 
   * @param xi 
   */
  void GetBu(const double &xi);

  /**
   * @Brief: Compute the Bw Vector
   * 
   * @param xi 
   */
  void GetBw(const double &xi);

  /**
   * @Brief: Compute the C Vector
   * 
   * @param xi 
   */
  void GetC(const double &xi);

  /**
   * @Brief: Transform Vector Variables into Global Coordinate System
   * 
   * @param aBar 
   * @param coord 
   * @return VectorXd 
   */
  VectorXd ToGlobalCoordinates(const VectorXd &aBar, const MatrixXd &coord);

  /**
   * @Brief:  Transform Matrix Variables into Global Coordinate System
   * 
   * @param ABar 
   * @param coord 
   * @return MatrixXd 
   */
  MatrixXd ToGlobalCoordinates(const MatrixXd &ABar, const MatrixXd &coord);

  /**
   * @Brief: Compute the Rotation Matrix
   * 
   * @param coord 
   * @return MatrixXd 
   */
  MatrixXd GetRotationMatrix(const MatrixXd &coord);

  /**
   * @Brief: Initialize Some Varibales
   * 
   */
  void Initialize();

private:
  double m_E = 0.;
  double m_A = 0.;
  double m_I = 0.;
  double m_G = 0.;
  double m_l0 = 0.;
  double m_EA = 0.;
  double m_EI = 0.;
  VectorXd m_Bu;
  VectorXd m_Bw;
  VectorXd m_C;

private:
  double m_epsl = 0.;
  double m_chi = 0.;
  double N = 0.;
  double M = 0.;
  double wght = 0.;
  double Jac = 0.;
  double tempDouble = 0.;
  VectorXd aBar;
  VectorXd a0;
  MatrixXd tempMat;
};
ReflectRegister(KirchhoffBeam, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif  // KIRCHHOFFBEAM_H