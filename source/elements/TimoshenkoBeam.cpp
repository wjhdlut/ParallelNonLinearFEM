#include <elements/TimoshenkoBeam.h>
#include <util/Tools.h>
#include <util/Transformations.h>

TimoshenkoBeam::TimoshenkoBeam(const std::vector<int> &elemNodes, const nlohmann::json &modelProps)
               : Element(elemNodes, modelProps)
{
  m_dofType = {"u", "v", "rz"};
  Tools::GetParameter(m_E, "E", modelProps);
  Tools::GetParameter(m_A, "A", modelProps);
  Tools::GetParameter(m_G, "G", modelProps);
  Tools::GetParameter(m_I, "I", modelProps);

  m_EA = m_E * m_A;
  m_EI = m_E * m_I;
  m_GA = m_G * m_A * 5. / 6.;

  m_intPoints = {std::make_pair<double, double>(-0.577, 1.0),
                 std::make_pair<double, double>( 0.577, 1.0)};
}

TimoshenkoBeam::~TimoshenkoBeam()
{
}

void TimoshenkoBeam::GetTangentStiffness(std::shared_ptr<ElementData> &elemDat)
{
  m_l0 = (elemDat->m_coords.row(2) - elemDat->m_coords.row(0)).norm();

  aBar = ToElemCoordinates(elemDat->m_state, elemDat->m_coords);

  for(auto intPoints : m_intPoints)
  {
    GetHt(intPoints.first);
    GetBu(intPoints.first);
    GetBw(intPoints.first);
    GetBt(intPoints.first);
    
    eps = bu.dot(aBar) + 0.5 * pow(bw.dot(aBar), 2.);
    gam = ht.dot(aBar) + bw.dot(aBar);
    chi = bt.dot(aBar);

    N = m_EA * eps;
    Q = m_GA * gam;
    M = m_EI * chi;

    wght = 0.5 * m_l0 * intPoints.second;

    // Compute internal force vector
    tempDouble = bw.dot(aBar);
    elemDat->m_fint += wght * N * bu + wght * (N * tempDouble + Q) * bw 
                     + wght * M * bt + wght * Q * ht;

    // Compute stiffness matrix
    elemDat->m_stiff += m_EA * wght * Math::VecCross(bu, bu) 
                      + m_EA * tempDouble * wght * (Math::VecCross(bu, bw) + Math::VecCross(bw, bu))
                      + wght * (m_EA * pow(tempDouble, 2) + m_GA + N) * Math::VecCross(bw, bw)
                      + wght * m_GA * (Math::VecCross(bw, ht) + Math::VecCross(ht, bw) + Math::VecCross(ht, ht))
                      + wght * m_EI * Math::VecCross(bt, bt);
  }

  elemDat->m_stiff(4, 4) = 1.;

  elemDat->m_fint = ToGlobalCoordinates(elemDat->m_fint, elemDat->m_coords);
  elemDat->m_stiff = ToGlobalCoordinates(elemDat->m_stiff, elemDat->m_coords);
}

void TimoshenkoBeam::GetHu(const double &xi)
{
  hu = VectorXd::Zero(9);

  hu(0) = 0.5 * (1. - xi);
  hu(3) = (1. - xi * xi);
  hu(6) = 0.5 * (1. + xi);
}

void TimoshenkoBeam::GetHw(const double &xi)
{
  hw = VectorXd::Zero(9);

  hw(1) = 0.5 * (1. - xi);
  hw(4) = (1. - xi * xi);
  hw(7) = 0.5 * (1. + xi);
}

void TimoshenkoBeam::GetHt(const double &xi)
{
  ht = VectorXd::Zero(9);

  ht(2) = 0.5 * (1. - xi);
  ht(5) = (1. - xi * xi);
  ht(8) = 0.5 * (1. + xi);
}

void TimoshenkoBeam::GetBu(const double &xi)
{
  bu = VectorXd::Zero(9);

  bu(0) = -1. / m_l0;
  bu(3) = -4. * xi / m_l0;
  bu(6) =  1. / m_l0;
}

void TimoshenkoBeam::GetBw(const double &xi)
{
  bw = VectorXd::Zero(9);

  bw(1) = -1. / m_l0;
  bw(7) =  1. / m_l0;

  bw(2) =  0.5 * xi;
  bw(8) = -0.5 * xi;
}

void TimoshenkoBeam::GetBt(const double &xi)
{
  bt = VectorXd::Zero(9);

  bt(2) = -1. / m_l0;
  bt(5) = -4. * xi / m_l0;
  bt(8) = 1. / m_l0;
}

VectorXd TimoshenkoBeam::ToElemCoordinates(const VectorXd&a, const MatrixXd &coords)
{
  MatrixXd R = GetRotationMatrix(coords);

  return R*a;
}

VectorXd TimoshenkoBeam::ToGlobalCoordinates(const VectorXd &aBar, const MatrixXd &coords)
{
  MatrixXd R = GetRotationMatrix(coords);

  return R.transpose() * aBar;
}

MatrixXd TimoshenkoBeam::ToGlobalCoordinates(const MatrixXd &ABar, const MatrixXd &coords)
{
  MatrixXd R = GetRotationMatrix(coords);

  return R.transpose() * ABar * R;
}


MatrixXd TimoshenkoBeam::GetRotationMatrix(const MatrixXd &coords)
{
  MatrixXd R = MatrixXd::Identity(9, 9);
  MatrixXd crd = MatrixXd::Zero(2, coords.cols());
  crd << coords.row(0), coords.row(2);

  MatrixXd tempMat = Transformations::GetRotationMatrix(crd);

  R(0, 0) = tempMat(0, 0), R(0, 1) = tempMat(0, 1);
  R(1, 0) = tempMat(1, 0), R(1, 1) = tempMat(1, 1);

  R(3, 3) = tempMat(0, 0), R(3, 4) = tempMat(0, 1);
  R(4, 3) = tempMat(1, 0), R(4, 4) = tempMat(1, 1);

  R(6, 6) = tempMat(0, 0), R(6, 7) = tempMat(0, 1);
  R(7, 6) = tempMat(1, 0), R(7, 7) = tempMat(1, 1);

  return R;
}