#ifndef FINITESTRAINCONTINUUM_H
#define FINITESTRAINCONTINUUM_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

/**
 * @Brief:  Total Lagrange Formulation
 * 
 */

class FiniteStrainContinuum : public Element
{
public:
  FiniteStrainContinuum(const std::string &elemShape,
                        const std::vector<int> &elemNodes,
                        const nlohmann::json &modelProps);
  
  ~FiniteStrainContinuum();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  void GetKinematics(const MatrixXd&dphi,
                     const VectorXd&elState);

  void GetBMatrix(const MatrixXd&dphi,
                  const MatrixXd&F);
  
  void Stress2Matrix(const VectorXd&stress);

  void GetBNLMatrix(const MatrixXd&dphi);

  void Initialize(const std::string &elemShape);

  // virtual void SetElemNodeOrdered();

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  
  double detJac = 0.;                                   // the Determinant of Jacobian Matrix
  /*jac = [pXpxi1 pXpxi2 pXpxi3,
           pYpxi1 pYpxi2 pYpxi3,
           pZpxi1 pZpxi2 pZpxi3];*/
  MatrixXd jac;                                         // the Jacobian Matrix
  MatrixXd invJac;                                      // the Inverse Jacobian Matrix

  /* pHpX = [pH1pX pH1pY pH1pZ,
             pH2pX pH2pY pH2pZ,
             ...
             pHnpX pHnpY pHnpZ]*/
  MatrixXd pHpX;                                        // the Derivative of Shape Function
  MatrixXd B;                                           // the Strain Matrix
  VectorXd sigma;                                       // the Stress Vector
  MatrixXd D;                                           // the Tangent Matrix
  MatrixXd Bnl;

  // stress matrix
  MatrixXd T; 
};

ReflectRegister(FiniteStrainContinuum, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUM_H