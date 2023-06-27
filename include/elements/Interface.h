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
  MatrixXd GetRotation(const MatrixXd &coords, const VectorXd &state);
  
  void GetBMatrix(const VectorXd &H, const MatrixXd &R);

  void GetKinematics(const VectorXd &elState);

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  MatrixXd B;                                 // the Strain Matrix
  VectorXd sigma;                             // the Stress Vector
  MatrixXd D;                                 // the Tangent Matrix
};

ReflectRegister(Interface, const std::vector<int> &, const nlohmann::json &)

#endif // INTERFACE_H