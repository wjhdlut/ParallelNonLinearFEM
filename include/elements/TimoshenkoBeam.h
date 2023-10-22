/**
 * @File Name:     TimoshenkoBeam.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-17
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef TIMOSHENKOBEAM_H
#define TIMOSHENKOBEAM_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class TimoshenkoBeam : public Element
{
public:
  /**
   * @Brief: Construct a new Timoshenko Beam object
   * 
   * @param elemShape 
   * @param elemNodes 
   * @param modelProps 
   */
  TimoshenkoBeam(const std::string &elemShape,
                 const std::vector<int> &elemNodes,
                 const nlohmann::json &modelProps);

  ~TimoshenkoBeam();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  void Initialize(const nlohmann::json &modelProps);

  void GetHu(const double &xi);

  void GetHw(const double &xi);

  void GetHt(const double &xi);

  void GetBu(const double &xi);

  void GetBw(const double &xi);

  void GetBt(const double &xi);

  VectorXd ToElemCoordinates(const VectorXd&a, const MatrixXd &coords);

  VectorXd ToGlobalCoordinates(const VectorXd &aBar, const MatrixXd &coords);

  MatrixXd ToGlobalCoordinates(const MatrixXd &ABar, const MatrixXd &coords);

  MatrixXd GetRotationMatrix(const MatrixXd &coords);

private:
  double m_E = 0.;
  double m_A = 0.;
  double m_I = 0.;
  double m_G = 0.;
  double m_l0 = 0.;
  double m_EA = 0.;
  double m_EI = 0.;
  double m_GA = 0.;
  std::vector<std::pair<double, double>> m_intPoints;


private:
  double tempDouble = 0.;
  double eps = 0., gam = 0., chi = 0.;
  double M = 0., Q = 0., N = 0.;
  double wght = 0.;
  VectorXd ht;
  VectorXd bu;
  VectorXd bw;
  VectorXd bt;
  VectorXd hu;
  VectorXd hw;
  VectorXd aBar;
};
ReflectRegister(TimoshenkoBeam, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif // TIMOSHENKOBEAM_H