
#ifndef SMALLSTRAINCONTINUUM_H
#define SMALLSTRAINCONTINUUM_H

#include <elements/Element.h>

class SmallStrainContinuum : public Element
{
public:
  SmallStrainContinuum(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  
  ~SmallStrainContinuum();

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;
};


#endif // SMALLSTRAINCONTINUUM_H