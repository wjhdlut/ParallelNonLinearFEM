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
  KirchhoffBeam(const std::vector<int> &elemNode, const nlohmann::json &modelProps);
  ~KirchhoffBeam();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  std::vector<double> TolemCoordinate(const std::vector<double> &a, const Matrix &coord);

  std::vector<double> GetBu(const double &xi);

  std::vector<double> GetBw(const double &xi);

  std::vector<double> GetC(const double &xi);

  std::vector<double> ToGlobalCoordinates(const std::vector<double> &aBar, const Matrix &coord);

  Matrix ToGlobalCoordinates(const Matrix &ABar, const Matrix &coord);

  Matrix GetRotationMatrix(const Matrix &coord);

private:
  double m_E = 0.;
  double m_A = 0.;
  double m_I = 0.;
  double m_G = 0.;
  double m_l0 = 0.;
  double m_EA = 0.;
  double m_EI = 0.;
  std::vector<double> m_weights = {5./9., 8./9., 5./9.};
  std::vector<double> m_intPoints = {-sqrt(3./5.), 0., sqrt(3./5.)};
  std::vector<double> m_Bu;
  std::vector<double> m_Bw;
  std::vector<double> m_C;

private:
  double m_epsl = 0.;
  double m_chi = 0.;
  double N = 0.;
  double M = 0.;
  double wght = 0.;
  double Jac = 0.;
  double tempDouble = 0.;
  std::vector<double> tempVec;
  Matrix tempMat;
};
ReflectRegister(KirchhoffBeam, const std::vector<int> &, const nlohmann::json &)

#endif  // KIRCHHOFFBEAM_H