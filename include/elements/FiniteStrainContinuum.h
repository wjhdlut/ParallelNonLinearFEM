#ifndef FINITESTRAINCONTINUUE_H
#define FINITESTRAINCONTINUUE_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class FiniteStrainContinuum : public Element
{
public:
  FiniteStrainContinuum(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  ~FiniteStrainContinuum();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

  virtual void GetMassMatrix(std::shared_ptr<ElementData> &elemDat) override;

private:
  void GetKinematics(const std::vector<std::vector<double>> &dphi,
                     const std::vector<double>&elState);

  void GetBMatrix(const std::vector<std::vector<double>>&dphi,
                  const std::vector<std::vector<double>>& F);
  
  void Stress2Matrix(const std::vector<double>&stress);

  void GetBNLMatrix(const Matrix&dphi);

  void GetNMatrix(const std::vector<double> &h);

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  
  double detJac = 0.;                                   // the Determinant of Jacobian Matrix
  /*jac = [pXpxi1 pXpxi2 pXpxi3,
             pYpxi1 pYpxi2 pYpxi3,
             pZpxi1 pZpxi2 pZpxi3];*/
  std::vector<std::vector<double>> jac;                  // the Jacobian Matrix
  std::vector<std::vector<double>> invJac;               // the Inverse Jacobian Matrix

  /* pHpX = [pH1pX1 pH1pX2 pH1pX3,
             pH2pX1 pH2pX2 pH2pX3,
             ...
             pHnpX1 pHnpX2 pHnpX3,]*/
  std::vector<std::vector<double>> pHpX;                 // the Derivative of Shape Function
  std::vector<std::vector<double>> B;                    // the Strain Matrix
  std::vector<double> sigma;                             // the Stress Vector
  std::vector<std::vector<double>> D;                    // the Tangent Matrix
  Matrix Bnl;

  // stress matrix
  std::vector<std::vector<double>> T; 

  Matrix N;
};

ReflectRegister(FiniteStrainContinuum, const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUE_H