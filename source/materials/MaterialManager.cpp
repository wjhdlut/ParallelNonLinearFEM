
#include <materials/MaterialManager.h>
#include <util/ObjectFactory.h>

MaterialManager::MaterialManager(const nlohmann::json &matProps)
{
  std::string matType = matProps.at("type");
  m_mat = ObjectFactory::CreateObject<BaseMaterial>(matType, matProps);
  iIter = -1;
}

MaterialManager::~MaterialManager()
{}

void MaterialManager::Reset()
{
  iIter = -1;
}