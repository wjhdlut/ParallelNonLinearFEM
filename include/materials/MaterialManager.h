#ifndef MATERIALMANAGER_H
#define MATERIALMANAGER_H

#include <nlohmann/json.hpp>
#include <materials/BaseMaterial.h>
#include <util/Kinematics.h>

class MaterialManager
{
public:
  MaterialManager(const nlohmann::json &matProps);
  ~MaterialManager();

  void Reset();

  std::vector<double> GetStress(const std::shared_ptr<Kinematics>&kin, int iSam = -1);

  inline std::vector<std::vector<double>> GetTangMatrix(){
    return m_mat->GetTangMatrix();
  }

private:
  int iIter = -1;
  std::shared_ptr<BaseMaterial> m_mat;
};

#endif // MATERIALMANAGER_H

