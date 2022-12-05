#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H

#include <nlohmann/json.hpp>
#include <materials/BaseMaterial.h>

class MaterialManager
{
public:
  MaterialManager(const nlohmann::json &matProps);
  ~MaterialManager();

  void Reset();

private:
  int iIter = -1;
  std::shared_ptr<BaseMaterial> m_mat;
};

#endif // MATERIALMANAGER_H

