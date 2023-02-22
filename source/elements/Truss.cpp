#include <elements/Truss.h>

Truss::Truss(const std::vector<int> &elemNode, const nlohmann::json &modelProps)
      : Element(elemNode, modelProps)
{

}

Truss::~Truss()
{}