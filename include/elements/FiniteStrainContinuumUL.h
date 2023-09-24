#ifndef FINITESTRAINCONTINUUMUL_H
#define FINITESTRAINCONTINUUMUL_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

/**
 * @Brief: Update Lagrange Formulation
 * 
 */

class FiniteStrainContinuumUL : public Element
{
public:
  FiniteStrainContinuumUL(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  ~FiniteStrainContinuumUL();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  void GetKinematics(const MatrixXd&dphi,
                     const VectorXd&elState);

  void GetBMatrix(const MatrixXd&dphi);
  
  void Stress2Matrix(const VectorXd&stress);

  void GetBNLMatrix(const MatrixXd&dphi);

  MatrixXd UpdateNodeCoords(std::shared_ptr<ElementData>&elemDat);

private:
  std::shared_ptr<Kinematics> kin = nullptr;
  
  /*jac = [pxpxi1 pxpxi2 pxpxi3,
           pypxi1 pypxi2 pypxi3,
           pzpxi1 pzpxi2 pzpxi3];*/
  MatrixXd jac;                                         // the Jacobian Matrix

  /*jac = [pXpxi1 pXpxi2 pXpxi3,
           pYpxi1 pYpxi2 pYpxi3,
           pZpxi1 pZpxi2 pZpxi3];*/
  MatrixXd Jac;                                         // the Jacobian Matrix
  
  MatrixXd invJac;                                      // the Inverse Jacobian Matrix

  /* pHpX = [pH1pX pH1pY pH1pZ,
             pH2pX pH2pY pH2pZ,
             ...
             pHnpX pHnpY pHnpZ]*/
  MatrixXd pHpX;

  /* pHpX = [pH1px pH1py pH1pz,
             pH2px pH2py pH2pz,
             ...
             pHnpx pHnpy pHnpz]*/
  MatrixXd pHpx;                                        // the Derivative of Shape Function
  MatrixXd B;                                           // the Strain Matrix
  VectorXd sigma;                                       // the Stress Vector
  MatrixXd D;                                           // the Tangent Matrix
  MatrixXd Bnl;
  MatrixXd T;                                           // stress matrix
  MatrixXd nodeCoord;
};

ReflectRegister(FiniteStrainContinuumUL, const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUMUL_H