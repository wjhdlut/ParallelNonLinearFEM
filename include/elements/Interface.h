#ifndef INTERFACE_H
#define INTERFACE_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class Interface : public Element
{
public:
  Interface(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  ~Interface();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  Matrix GetRotation(const Matrix &coords, const std::vector<double> &state);
  
  void GetBMatrix(const std::vector<double> &H, const Matrix &R);

  void GetKinematics(const std::vector<double> &elState);

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  std::vector<std::vector<double>> B;                    // the Strain Matrix
  std::vector<double> sigma;                             // the Stress Vector
  std::vector<std::vector<double>> D;                    // the Tangent Matrix
};

ReflectRegister(Interface, const std::vector<int> &, const nlohmann::json &)

#endif // INTERFACE_H