
#include <elements/FiniteStrainContinuum.h>

FiniteStrainContinuum::FiniteStrainContinuum(const std::vector<int> &elemNodes,
                                             const nlohmann::json &modelProps)
                      : Element(elemNodes, modelProps)
{
  m_dofType = {"u", "v"};
}

FiniteStrainContinuum::~FiniteStrainContinuum()
{}

void FiniteStrainContinuum::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{
  
}

