

#include <elements/SmallStrainContinuum.h>

SmallStrainContinuum::SmallStrainContinuum(const std::vector<int> &elemNodes,
                                           const nlohmann::json &modelProps)
                                          : Element(elemNodes, modelProps)
{
}

SmallStrainContinuum::~SmallStrainContinuum()
{}

void SmallStrainContinuum::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{

}