#include <elements/Spring.h>
#include <util/Transformations.h>
#include <util/Tools.h>

Spring::Spring(const std::vector<int> &elemNode, const nlohmann::json &modelProps)
       : Element(elemNode, modelProps)
{
  Tools::GetParameter(m_k, "k", m_props);
}

Spring::~Spring()
{}

void Spring::GetTangentStiffness(std::shared_ptr<ElementData> &elemDat)
{
  // Compute the current state vector
  std::vector<double> a = Transformations::ToElementCoordinates(elemDat->m_state, elemDat->m_coords);
  std::vector<double> Da = Transformations::ToElementCoordinates(elemDat->m_Dstate, elemDat->m_coords);

  // Compute the elongation of the spring
  double elong = a[2] - a[0];

  // Compute the force in the spring
  double Fs = elong * m_k;

  // Compute the element internal force vector in the element coordinate system
  std::vector<double> elemFint = {-Fs, 0., Fs, 0};

  // Determine the element tangent stiffness in the element coordinate system
  std::vector<double> tempVec(4, 0.);
  std::vector<std::vector<double>> elemK(4, tempVec);
  elemK[0][0] = m_k;
  elemK[1][1] = m_k;

  elemK[0][2] = -m_k;
  elemK[1][3] = -m_k;

  elemK[2][0] = -m_k;
  elemK[3][1] = -m_k;

  elemK[2][2] = m_k;
  elemK[3][3] = m_k;

  elemDat->m_stiff = Transformations::ToGlobalCoordinates(elemK, elemDat->m_coords);
  elemDat->m_fint = Transformations::ToGlobalCoordinates(elemFint, elemDat->m_coords);
}