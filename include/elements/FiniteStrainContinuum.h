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
  void GetKinematics(const MatrixXd&dphi,
                     const VectorXd&elState);

  void GetBMatrix(const MatrixXd&dphi,
                  const MatrixXd&F);
  
  void Stress2Matrix(const VectorXd&stress);

  void GetBNLMatrix(const MatrixXd&dphi);

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  
  double detJac = 0.;                                   // the Determinant of Jacobian Matrix
  /*jac = [pXpxi1 pXpxi2 pXpxi3,
           pYpxi1 pYpxi2 pYpxi3,
           pZpxi1 pZpxi2 pZpxi3];*/
  MatrixXd jac;                                         // the Jacobian Matrix
  MatrixXd invJac;                                      // the Inverse Jacobian Matrix

  /* pHpX = [pH1pX1 pH1pX2 pH1pX3,
             pH2pX1 pH2pX2 pH2pX3,
             ...
             pHnpX1 pHnpX2 pHnpX3,]*/
  MatrixXd pHpX;                                        // the Derivative of Shape Function
  MatrixXd B;                                           // the Strain Matrix
  VectorXd sigma;                                       // the Stress Vector
  MatrixXd D;                                           // the Tangent Matrix
  MatrixXd Bnl;

  // stress matrix
  MatrixXd T; 
};

ReflectRegister(FiniteStrainContinuum, const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUE_H