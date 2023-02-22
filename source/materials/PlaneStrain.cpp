
#include <materials/PlaneStrain.h>


PlaneStrain::PlaneStrain(const nlohmann::json &props) : BaseMaterial(props)
{
  m_E = props.at("E"); m_nu = props.at("nu");

  ComputeDMatrix();
}

PlaneStrain::~PlaneStrain()
{}

void PlaneStrain::ComputeDMatrix()
{
  std::vector<double> temp(3, 0.);
  m_D.resize(3, temp);

  m_D[0][0] = m_E*(1.-m_nu)/((1+m_nu)*(1.-2.*m_nu));
  m_D[0][1] = m_D[0][0]*m_nu/(1-m_nu);
  m_D[1][0] = m_D[0][1];
  m_D[1][1] = m_D[0][0];
  m_D[2][2] = m_D[0][0]*0.5*(1.-2.*m_nu)/(1.-m_nu);
}
