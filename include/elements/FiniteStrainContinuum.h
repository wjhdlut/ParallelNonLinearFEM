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

private:
  std::shared_ptr<Kinematics> GetKinematics(const std::vector<std::vector<double>> &dphi,
                                            const std::vector<double>&elState);

  std::vector<std::vector<double>> GetBMatrix(const std::vector<std::vector<double>>&dphi,
                                              const std::vector<std::vector<double>>& F);
  
  std::vector<std::vector<double>> Stress2Matrix(const std::vector<double>&stress);

  Matrix GetBNLMatrix(const Matrix&dphi);

private:
  int order = 0;
  double detJac = 0.;
  /*jac = [pXpxi1 pXpxi2 pXpxi3,
             pYpxi1 pYpxi2 pYpxi3,
             pZpxi1 pZpxi2 pZpxi3];*/
  std::vector<std::vector<double>> jac;
  std::vector<std::vector<double>> invJac;
  
  /* pHpX = [pH1pX1 pH1pX2 pH1pX3,
             pH2pX1 pH2pX2 pH2pX3,
             ...
             pHnpX1 pHnpX2 pHnpX3,]*/
  std::vector<std::vector<double>> pHpX;

  // Strain Matrix
  std::vector<std::vector<double>> B;
  Matrix Bnl;

  // stress vector
  std::vector<double> sigma;

  // stress matrix
  std::vector<std::vector<double>> T;

  // tangent Matrix
  std::vector<std::vector<double>> D;
  std::shared_ptr<Kinematics> kin = nullptr;

  Matrix tempMatrix;
  double weighti = 0.;
  std::string method = "Gauss";

  std::vector<std::vector<double>> xi;
  std::vector<double> weight;

  std::vector<double> tempVec;

  Matrix outputData;
};

ReflectRegister(FiniteStrainContinuum, const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUE_H