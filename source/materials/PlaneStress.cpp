#include <materials/PlaneStress.h>
#include <util/Math.h>

PlaneStress::PlaneStress(const nlohmann::json &props) : BaseMaterial(props)
{
  SetHistoryParameter("kappa", 0.);
  CommitHistory();
  
  if(m_props.at("E").is_string()){
    std::string E = m_props.at("E");
    m_E = std::stod(E);
  }
  else{
    m_E = m_props.at("E");
  }

  if(m_props.at("nu").is_string()){
    std::string nu = m_props.at("nu");
    m_nu = std::stod(nu);
  }
  else{
    m_nu = m_props.at("nu");
  }
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