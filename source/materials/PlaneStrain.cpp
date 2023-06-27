
#include <materials/PlaneStrain.h>


PlaneStrain::PlaneStrain(const nlohmann::json &props) : BaseMaterial(props)
{
  m_E = SetMaterialParamter("E");
  m_nu = SetMaterialParamter("nu");
  m_rho = SetMaterialParamter("rho");
  ComputeDMatrix();
}

PlaneStrain::~PlaneStrain()
{}

void PlaneStrain::ComputeDMatrix()
{
  m_D = MatrixXd::Zero(3, 3);

  m_D(0, 0) = m_E*(1.-m_nu)/((1+m_nu)*(1.-2.*m_nu));
  m_D(0, 1) = m_D(0, 0)*m_nu/(1-m_nu);
  m_D(1, 0) = m_D(0, 1);
  m_D(1, 1) = m_D(0, 0);
  m_D(2, 2) = m_D(0, 0)*0.5*(1.-2.*m_nu)/(1.-m_nu);
}
