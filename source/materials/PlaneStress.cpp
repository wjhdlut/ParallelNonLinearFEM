#include <materials/PlaneStress.h>

PlaneStress::PlaneStress(const nlohmann::json &props) : BaseMaterial(props)
{
  SetHistoryParameter("kappa", 0.);
  CommitHistory();
  
  m_E = m_props.at("E"); m_nu = m_props.at("nu");
  ComputeDMatrix();
}

PlaneStress::~PlaneStress()
{}

void PlaneStress::ComputeDMatrix()
{
  m_D[0][0] = m_E/(1. - m_nu * m_nu);
  m_D[0][1] = m_D[0][0]*m_nu;
  m_D[1][0] = m_D[0][1];
  m_D[1][1] = m_D[0][0];
  m_D[2][2] = m_E/2./(1. + m_nu);
}
