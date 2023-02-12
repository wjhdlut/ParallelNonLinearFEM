#include <materials/PlaneStress.h>
#include <util/Math.h>

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
  std::vector<double> temp(3, 0.);
  m_D.resize(3, temp);
  m_D[0][0] = m_E/(1. - m_nu * m_nu);
  m_D[0][1] = m_D[0][0]*m_nu;
  m_D[1][0] = m_D[0][1];
  m_D[1][1] = m_D[0][0];
  m_D[2][2] = m_E/2./(1. + m_nu);
}

std::vector<double> PlaneStress::GetStress(const std::shared_ptr<Kinematics>&kin){
  return Math::MatraixAMultVecB(m_D, kin->strain);
}