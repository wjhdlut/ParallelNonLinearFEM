
#include <materials/LinearElasticity.h>

LinearElasticity::LinearElasticity(const nlohmann::json &props) : BaseMaterial(props)
{
  m_E = SetMaterialParamter("E");
  m_nu = SetMaterialParamter("nu");
  ComputeDMatrix();
}

LinearElasticity::~LinearElasticity()
{}

void LinearElasticity::ComputeDMatrix()
{
  double fac = 1.0 / (2.0 * m_nu * m_nu + m_nu - 1.0 );
  
  m_D = MatrixXd::Zero(6, 6);

  m_D(0, 0) = fac * m_E * ( m_nu - 1.0 );
  m_D(0, 1) = -1.0 * fac * m_E * m_nu;
  m_D(0, 2) = m_D(0, 1);
  m_D(1, 0) = m_D(0, 1);
  m_D(1, 1) = m_D(0, 0);
  m_D(1, 2) = m_D(0, 1);
  m_D(2, 0) = m_D(0, 1);
  m_D(2, 1) = m_D(0, 1);
  m_D(2, 2) = m_D(0, 0);
  m_D(3, 3) = m_E / ( 2.0 + 2.0 * m_nu );
  m_D(4, 4) = m_D(3, 3);
  m_D(5, 5) = m_D(3, 3);
}