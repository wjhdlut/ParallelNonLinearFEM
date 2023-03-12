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
  std::shared_ptr<Kinematics> GetKinematics(const std::vector<std::vector<double>> &dphi,
                                            const std::vector<double>&elState);

  std::vector<std::vector<double>> GetBMatrix(const std::vector<std::vector<double>>&dphi,
                                              const std::vector<std::vector<double>>& F);
  
  std::vector<std::vector<double>> Stress2Matrix(const std::vector<double>&stress);

  Matrix GetBNLMatrix(const Matrix&dphi);

  Matrix GetNMatrix(const std::vector<double> &h);

private:
  Matrix Bnl;

  // stress matrix
  std::vector<std::vector<double>> T; 

  // std::vector<std::vector<double>> xi;
  // std::vector<double> weight;
};

ReflectRegister(FiniteStrainContinuum, const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUE_H