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
};

ReflectRegister(FiniteStrainContinuum, const std::vector<int> &, const nlohmann::json &)

#endif // FINITESTRAINCONTINUUE_H