#ifndef TIMOSHENKOBEAM_H
#define TIMOSHENKOBEAM_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class TimoshenkoBeam : public Element
{
public:
  TimoshenkoBeam(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  ~TimoshenkoBeam();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  std::vector<double> GetHu(const double &xi);

  std::vector<double> GetHw(const double &xi);

  std::vector<double> GetHt(const double &xi);

  std::vector<double> GetBu(const double &xi);

  std::vector<double> GetBw(const double &xi);

  std::vector<double> GetBt(const double &xi);

  std::vector<double> ToElemCoordinates(const std::vector<double>&a, const Matrix &coords);

  std::vector<double> ToGlobalCoordinates(const std::vector<double> &aBar, const Matrix &coords);

  Matrix ToGlobalCoordinates(const Matrix &ABar, const Matrix &coords);

  Matrix GetRotationMatrix(const Matrix &coords);

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
  std::vector<double> ht;
  std::vector<double> bu;
  std::vector<double> bw;
  std::vector<double> bt;
};
ReflectRegister(TimoshenkoBeam, const std::vector<int> &, const nlohmann::json &)

#endif // TIMOSHENKOBEAM_H