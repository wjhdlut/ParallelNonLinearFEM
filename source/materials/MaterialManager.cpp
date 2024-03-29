
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

VectorXd MaterialManager::GetStress(const std::shared_ptr<Kinematics>&kin,
                                               const VectorXd &increDisp,
                                               const MatrixXd &dphi, int iSam)
{
  if(-1 == iSam){
    iIter += 1;
    iSam = iIter;
  }

  m_mat->SetIter(iSam);

  return m_mat->GetStress(kin, increDisp, dphi);
}