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
  
  weight = VectorXd::Zero(3);
  weight << 5./9., 8./9., 5./9.;
  
  xi = MatrixXd::Zero(3, 1);
  xi << -sqrt(3./5.), 0., sqrt(3./5.);
}

KirchhoffBeam::~KirchhoffBeam()
{}

void KirchhoffBeam::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{
  a0 = elemDat->m_coords.row(2) - elemDat->m_coords.row(0);
  m_l0 = a0.norm();
  Jac =  0.5 * m_l0;

  aBar = TolemCoordinate(elemDat->m_state, elemDat->m_coords);
  
  for(int i = 0; i < xi.rows(); i++)
  {
    GetBu(xi(i, 0));
    GetBw(xi(i, 0));
    GetC(xi(i, 0));
    
    m_epsl = m_Bu.dot(aBar) + 0.5 * pow(m_Bw.dot(aBar), 2.);
    m_chi  = m_C.dot(aBar);

    N = m_EA * m_epsl;
    M = m_EI * m_chi;
    wght = Jac * weight(i);
    tempDouble = m_Bw.dot(aBar);

    // Compute internale force vector
    elemDat->m_fint += N * wght * m_Bu + N * wght * tempDouble * m_Bw + wght * M * m_C;

    // Compute stiffness matrix
    elemDat->m_stiff += m_EA * wght * Math::VecCross(m_Bu, m_Bu)
                      + m_EA * tempDouble * wght * (Math::VecCross(m_Bu, m_Bw) + Math::VecCross(m_Bw, m_Bu))
                      + m_EI * wght * Math::VecCross(m_C, m_C)
                      + wght * (m_EA * pow(tempDouble, 2.) + N) * Math::VecCross(m_Bw, m_Bw);
  }
  elemDat->m_stiff(4, 4) = 1., elemDat->m_stiff(5, 5) = 1.;

  // Transform from element local coordinate system to global coordinate system
  elemDat->m_fint = ToGlobalCoordinates(elemDat->m_fint, elemDat->m_coords);
  elemDat->m_stiff = ToGlobalCoordinates(elemDat->m_stiff, elemDat->m_coords);
}

VectorXd KirchhoffBeam::TolemCoordinate(const VectorXd &a, const MatrixXd &coord)
{
  MatrixXd R = GetRotationMatrix(coord);

  return R * a;
}

void KirchhoffBeam::GetBu(const double &xi)
{
  int number = pow(m_dofType.size(), 2);
  m_Bu = VectorXd::Zero(number);
  
  m_Bu(0) = -1./m_l0;
  m_Bu(3) = -4.*xi/m_l0;
  m_Bu(6) =  1./m_l0;
}

void KirchhoffBeam::GetBw(const double &xi)
{
  int number = pow(m_dofType.size(), 2);
  m_Bw = VectorXd::Zero(number);
  
  m_Bw(1) = 1.5/m_l0*(xi*xi-1.0);
  m_Bw(2) = 0.25*(3*xi*xi-2.0*xi-1.0);
  m_Bw(7) = -1.5/m_l0*(xi*xi-1.0);
  m_Bw(8) = 0.25*(3*xi*xi+2.0*xi-1.0);
}

void KirchhoffBeam::GetC(const double &xi)
{
  int number = pow(m_dofType.size(), 2);
  m_C = VectorXd::Zero(number);
  
  m_C(1) = 6.0*xi/m_l0/m_l0;
  m_C(2) = (3.0*xi-1.0)/m_l0;
  m_C(7) = -6.0*xi/m_l0/m_l0;
  m_C(8) = (3.0*xi+1.0)/m_l0;
}

VectorXd KirchhoffBeam::ToGlobalCoordinates(const VectorXd &aBar, const MatrixXd &coord)
{
  MatrixXd R = GetRotationMatrix(coord);

  return R.transpose() * aBar;
}

MatrixXd KirchhoffBeam::ToGlobalCoordinates(const MatrixXd &ABar, const MatrixXd &coord)
{
  MatrixXd R = GetRotationMatrix(coord);

  return R.transpose() * ABar * R;
}

MatrixXd KirchhoffBeam::GetRotationMatrix(const MatrixXd &coord)
{
  MatrixXd R = MatrixXd::Identity(9, 9);

  MatrixXd crd = MatrixXd::Zero(2, coord.cols());
  crd << coord.row(0), coord.row(2);

  MatrixXd tempR = Transformations::GetRotationMatrix(crd);

  R(0, 0) = tempR(0, 0), R(0, 1) = tempR(0, 1);
  R(1, 0) = tempR(1, 0), R(1, 1) = tempR(1, 1);

  R(6, 6) = tempR(0, 0), R(6, 7) = tempR(0, 1);
  R(7, 6) = tempR(1, 0), R(7, 7) = tempR(1, 1);

  return R;
}