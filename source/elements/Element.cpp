
#include <elements/Element.h>

Element::Element(const std::vector<int> &elemNodes, const nlohmann::json &modelProps)
        : m_nodes(elemNodes)
{

  for(auto iter = modelProps.begin(); iter != modelProps.end(); iter++)
  {
    if("material" == iter.key())
    {
      const nlohmann::json &matProps = iter.value();
      m_mat = std::make_shared<MaterialManager>(matProps);
    }
    else{
      m_props[iter.key()] = iter.value();
    }
  }
}

Element::~Element()
{
}
