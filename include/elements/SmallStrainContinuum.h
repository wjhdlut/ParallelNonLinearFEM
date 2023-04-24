
#ifndef SMALLSTRAINCONTINUUM_H
#define SMALLSTRAINCONTINUUM_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>

class SmallStrainContinuum : public Element
{
public:
  SmallStrainContinuum(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  
  ~SmallStrainContinuum();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  void ComputeElemTimeStep(const std::shared_ptr<ElementShapeFunctions> &res,
                           const std::shared_ptr<ElementData> &elemDat);
  
  void GetKinematics(const std::vector<double> &elState);

  void GetBMatrix(const Matrix &dphi);

  void HourGlassTech(std::shared_ptr<ElementData>&elemDat,
                     const std::shared_ptr<ElementShapeFunctions> &res);

private:
  double m_vol = 0;
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
};

ReflectRegister(SmallStrainContinuum, const std::vector<int> &, const nlohmann::json &)

#endif // SMALLSTRAINCONTINUUM_H