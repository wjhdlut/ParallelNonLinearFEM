#include <elements/Spring.h>
#include <util/Transformations.h>
#include <util/Tools.h>
#include <iostream>

Spring::Spring(const std::string &elemShape,
               const std::vector<int> &elemNode,
               const nlohmann::json &modelProps)
       : Element(elemNode, modelProps)
{
  Tools::GetParameter(m_k, "k", m_props);
}

Spring::~Spring()
{}

void Spring::GetTangentStiffness(std::shared_ptr<ElementData> &elemDat)
{
  // Compute the current state vector
  a = Transformations::ToElementCoordinates(elemDat->m_state, elemDat->m_coords);
  Da = Transformations::ToElementCoordinates(elemDat->m_Dstate, elemDat->m_coords);

  // Compute the elongation of the spring
  elong = a[2] - a[0];

  // Compute the force in the spring
  Fs = elong * m_k;

  // Compute the element internal force vector in the element coordinate system
  elemDat->m_fint << -Fs, 0., Fs, 0;

  // Determine the element tangent stiffness in the element coordinate system
  elemDat->m_stiff(0, 0) = m_k;
  elemDat->m_stiff(1, 1) = m_k;

  elemDat->m_stiff(0, 2) = -m_k;
  elemDat->m_stiff(1, 3) = -m_k;

  elemDat->m_stiff(2, 0) = -m_k;
  elemDat->m_stiff(3, 1) = -m_k;

  elemDat->m_stiff(2, 2) = m_k;
  elemDat->m_stiff(3, 3) = m_k;

  elemDat->m_stiff = Transformations::ToGlobalCoordinates(elemDat->m_stiff, elemDat->m_coords);
  elemDat->m_fint = Transformations::ToGlobalCoordinates(elemDat->m_fint, elemDat->m_coords);
}