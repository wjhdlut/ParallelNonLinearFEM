/**
 * @File Name:     KirchhoffBeam.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-03-13
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <elements/KirchhoffBeam.h>
#include <util/Transformations.h>
#include <util/Tools.h>

KirchhoffBeam::KirchhoffBeam(const std::vector<int>&elemNode,
                             const nlohmann::json &modelProps)
              : Element(elemNode, modelProps)
{
  m_dofType = {"u", "v", "rz"};

  Tools::GetParameter(m_E, "E", m_props);
  Tools::GetParameter(m_A, "A", m_props);
  Tools::GetParameter(m_I, "I", m_props);
  Tools::GetParameter(m_G, "G", m_props);

  m_EI = m_E * m_I;
  m_EA = m_E * m_A;
}

KirchhoffBeam::~KirchhoffBeam()
{}

void KirchhoffBeam::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{
  std::vector<double> a0 = Math::VecAdd(-1, elemDat->m_coords[2], elemDat->m_coords[0]);
  m_l0 = Math::VecNorm(a0);
  Jac =  0.5 * m_l0;

  std::vector<double> aBar = TolemCoordinate(elemDat->m_state, elemDat->m_coords);
  
  int count = 0;
  for(auto xi : m_intPoints)
  {
    m_Bu = GetBu(xi);
    m_Bw = GetBw(xi);
    m_C = GetC(xi);

    m_epsl = Math::VecDot(m_Bu, aBar) + 0.5 * pow(Math::VecDot(m_Bw, aBar), 2.);
    m_chi = Math::VecDot(m_C, aBar);

    N = m_EA * m_epsl;
    M = m_EI * m_chi;

    wght = Jac * m_weights[count];

    tempDouble = Math::VecDot(m_Bw, aBar);
    // Compute internale force vector
    elemDat->m_fint = Math::VecAdd(N*wght, elemDat->m_fint, m_Bu);
    tempVec = Math::VecAdd(N/M*tempDouble, m_C, m_Bw);
    elemDat->m_fint = Math::VecAdd(wght*M, elemDat->m_fint, tempVec);

    // Compute stiffness matrix
    elemDat->m_stiff = Math::MatrixAdd(m_EA*wght, elemDat->m_stiff, Math::VecOuter(m_Bu, m_Bu));
    elemDat->m_stiff = Math::MatrixAdd(m_EA*tempDouble*wght, elemDat->m_stiff, Math::VecOuter(m_Bu, m_Bw));
    elemDat->m_stiff = Math::MatrixAdd(m_EA*tempDouble*wght, elemDat->m_stiff, Math::VecOuter(m_Bw, m_Bu));

    tempMat = Math::VecOuter(m_C, m_C);
    elemDat->m_stiff = Math::MatrixAdd(m_EI * wght, elemDat->m_stiff, tempMat);

    tempMat = Math::VecOuter(m_Bw, m_Bw);
    elemDat->m_stiff = Math::MatrixAdd(m_EA * pow(tempDouble, 2.) * wght, elemDat->m_stiff, tempMat);

    tempMat = Math::VecOuter(m_Bw, m_Bw);
    elemDat->m_stiff = Math::MatrixAdd(N * wght, elemDat->m_stiff, tempMat);

    elemDat->m_fint;

    count += 1;
  }
  elemDat->m_stiff[4][4] = 1., elemDat->m_stiff[5][5] = 1.;

  // Transform from element local coordinate system to global coordinate system
  elemDat->m_fint = ToGlobalCoordinates(elemDat->m_fint, elemDat->m_coords);
  elemDat->m_stiff = ToGlobalCoordinates(elemDat->m_stiff, elemDat->m_coords);
}

std::vector<double> KirchhoffBeam::TolemCoordinate(const std::vector<double> &a,
                                                   const Matrix &coord)
{
  Matrix R = GetRotationMatrix(coord);

  return Math::MatrixAMultVecB(R, a);
}

std::vector<double> KirchhoffBeam::GetBu(const double &xi)
{
  int number = pow(m_dofType.size(), 2);
  std::vector<double> Bu(number, 0);
  
  Bu[0] = -1./m_l0;
  Bu[3] = -4.*xi/m_l0;
  Bu[6] =  1./m_l0;

  return Bu;
}

std::vector<double> KirchhoffBeam::GetBw(const double &xi)
{
  int number = pow(m_dofType.size(), 2);
  std::vector<double> Bw(number, 0);
  
  Bw[1] = 1.5/m_l0*(xi*xi-1.0);
  Bw[2] = 0.25*(3*xi*xi-2.0*xi-1.0);
  Bw[7] = -1.5/m_l0*(xi*xi-1.0);
  Bw[8] = 0.25*(3*xi*xi+2.0*xi-1.0);
  return Bw;
}

std::vector<double> KirchhoffBeam::GetC(const double &xi)
{
  int number = pow(m_dofType.size(), 2);
  std::vector<double> C(number, 0);
  
  C[1] = 6.0*xi/m_l0/m_l0;
  C[2] = (3.0*xi-1.0)/m_l0;
  C[7] = -6.0*xi/m_l0/m_l0;
  C[8] = (3.0*xi+1.0)/m_l0;

  return C;
}

std::vector<double> KirchhoffBeam::ToGlobalCoordinates(const std::vector<double> &aBar, const Matrix &coord)
{
  Matrix R = GetRotationMatrix(coord);

  return Math::MatrixATransMultVecB(R, aBar);
}

Matrix KirchhoffBeam::ToGlobalCoordinates(const Matrix &ABar, const Matrix &coord)
{
  Matrix R = GetRotationMatrix(coord);

  tempMat = Math::MatrixAMultB(ABar, R);

  return Math::MatrixATransMultB(R, tempMat);
}

Matrix KirchhoffBeam::GetRotationMatrix(const Matrix &coord)
{
  Matrix R = Math::MatrixEye(9);

  Matrix crd = {coord[0], coord[2]};

  Matrix tempR = Transformations::GetRotationMatrix(crd);

  R[0][0] = tempR[0][0], R[0][1] = tempR[0][1], R[0][2] = tempR[0][2];
  R[1][0] = tempR[1][0], R[1][1] = tempR[1][1], R[1][2] = tempR[1][2];
  R[2][0] = tempR[2][0], R[2][1] = tempR[2][1], R[2][2] = tempR[2][2];

  R[6][6] = tempR[0][0], R[6][7] = tempR[0][1], R[6][8] = tempR[0][2];
  R[7][6] = tempR[1][0], R[7][7] = tempR[1][1], R[7][8] = tempR[1][2];
  R[8][6] = tempR[2][0], R[8][7] = tempR[2][1], R[8][8] = tempR[2][2];

  return R;
}