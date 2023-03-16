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

  m_intPoints = {std::make_pair<double, double>(-0.5774, 1.0),
                 std::make_pair<double, double>(0.5774, 1.0)};
}

TimoshenkoBeam::~TimoshenkoBeam()
{
}

void TimoshenkoBeam::GetTangentStiffness(std::shared_ptr<ElementData> &elemDat)
{
  m_l0 = Math::VecNorm(Math::VecAdd(-1., elemDat->m_coords[2], elemDat->m_coords[0]));

  std::vector<double> aBar = ToElemCoordinates(elemDat->m_state, elemDat->m_coords);

  for(auto intPoints : m_intPoints)
  {
    ht = GetHt(intPoints.first);
    bu = GetBu(intPoints.first);
    bw = GetBw(intPoints.first);
    bt = GetBt(intPoints.first);

    eps = Math::VecDot(bu, aBar) + 0.5 * pow(Math::VecDot(bw, aBar), 2.);
    gam = Math::VecDot(ht, aBar) + Math::VecDot(bw, aBar);
    chi = Math::VecDot(bt, aBar);

    N = m_EA * eps;
    Q = m_GA * gam;
    N = m_EI * chi;

    wght = 0.5 * m_l0 * intPoints.second;

    // Compute internal force vector
    tempDouble = Math::VecDot(bw, aBar);
    elemDat->m_fint = Math::VecAdd(wght*N, elemDat->m_fint, bu);
    elemDat->m_fint = Math::VecAdd(wght*(N*tempDouble + Q), elemDat->m_fint, bw);
    elemDat->m_fint = Math::VecAdd(wght*M, elemDat->m_fint, bt);
    elemDat->m_fint = Math::VecAdd(wght*Q, elemDat->m_fint, ht);

    // Compute stiffness matrix
    elemDat->m_stiff = Math::MatrixAdd(m_EA * wght, elemDat->m_stiff, Math::VecOuter(bu, bu));
    elemDat->m_stiff = Math::MatrixAdd(m_EA * tempDouble * wght, elemDat->m_stiff, Math::VecOuter(bu, bw));
    elemDat->m_stiff = Math::MatrixAdd(m_EA * tempDouble * wght, elemDat->m_stiff, Math::VecOuter(bw, bu));
    elemDat->m_stiff = Math::MatrixAdd(wght * (m_EA * pow(tempDouble, 2) + m_GA + N),
                                       elemDat->m_stiff, Math::VecOuter(bw, bw));
    elemDat->m_stiff = Math::MatrixAdd(wght * m_GA, elemDat->m_stiff, Math::VecOuter(bw, ht));
    elemDat->m_stiff = Math::MatrixAdd(wght * m_GA, elemDat->m_stiff, Math::VecOuter(ht, bw));
    elemDat->m_stiff = Math::MatrixAdd(wght * m_EI, elemDat->m_stiff, Math::VecOuter(bw, ht));
    elemDat->m_stiff = Math::MatrixAdd(wght * m_GA, elemDat->m_stiff, Math::VecOuter(ht, ht));
  }

  elemDat->m_stiff[4][4] = 1.;

  elemDat->m_fint = ToGlobalCoordinates(elemDat->m_fint, elemDat->m_coords);
  elemDat->m_stiff = ToGlobalCoordinates(elemDat->m_stiff, elemDat->m_coords);
}

std::vector<double> TimoshenkoBeam::GetHu(const double &xi)
{
  std::vector<double> Hu(9, 0.);

  Hu[0] = 0.5 * (1. - xi);
  Hu[3] = (1. - xi * xi);
  Hu[6] = 0.5 * (1. + xi);

  return Hu;
}

std::vector<double> TimoshenkoBeam::GetHw(const double &xi)
{
  std::vector<double> Hw(9, 0.);

  Hw[1] = 0.5 * (1. - xi);
  Hw[4] = (1. - xi * xi);
  Hw[7] = 0.5 * (1. + xi);

  return Hw;
}

std::vector<double> TimoshenkoBeam::GetHt(const double &xi)
{
  std::vector<double> Ht(9, 0.);

  Ht[2] = 0.5 * (1. - xi);
  Ht[5] = (1. - xi * xi);
  Ht[8] = 0.5 * (1. + xi);

  return Ht;
}

std::vector<double> TimoshenkoBeam::GetBu(const double &xi)
{
  std::vector<double> Bu(9, 0.);

  Bu[0] = -1. / m_l0;
  Bu[3] = -4. * xi / m_l0;
  Bu[6] =  1. / m_l0;

  return Bu;
}

std::vector<double> TimoshenkoBeam::GetBw(const double &xi)
{
  std::vector<double> Bw(9, 0.);

  Bw[1] = -1. / m_l0;
  Bw[7] =  1. / m_l0;

  Bw[2] =  0.5 * xi;
  Bw[8] = -0.5 * xi;
 
  return Bw;
}

std::vector<double> TimoshenkoBeam::GetBt(const double &xi)
{
  std::vector<double> Bt(9, 0.);

  Bt[2] = -1. / m_l0;
  Bt[5] = -4. * xi / m_l0;
  Bt[8] = 1. / m_l0;

  return Bt;
}

std::vector<double> TimoshenkoBeam::ToElemCoordinates(const std::vector<double>&a, 
                                                      const Matrix &coords)
{
  Matrix R = GetRotationMatrix(coords);

  return Math::MatrixAMultVecB(R, a);
}

std::vector<double> TimoshenkoBeam::ToGlobalCoordinates(const std::vector<double> &aBar, 
                                                        const Matrix &coords)
{
  Matrix R = GetRotationMatrix(coords);

  return Math::MatrixATransMultVecB(R, aBar);
}

Matrix TimoshenkoBeam::ToGlobalCoordinates(const Matrix &ABar, const Matrix &coords)
{
  Matrix R = GetRotationMatrix(coords);

  return Math::MatrixATransMultB(R, Math::MatrixAMultB(ABar, R));
}


Matrix TimoshenkoBeam::GetRotationMatrix(const Matrix &coords)
{
  Matrix R = Math::MatrixEye(9);
  Matrix crd = {coords[0], coords[2]};

  Matrix tempMat = Transformations::GetRotationMatrix(crd);

  R[3][3] = tempMat[0][0], R[3][4] = tempMat[0][1], R[3][5] = tempMat[0][2];
  R[4][3] = tempMat[1][0], R[4][4] = tempMat[1][1], R[4][5] = tempMat[1][2];
  R[5][3] = tempMat[2][0], R[3][4] = tempMat[2][1], R[3][5] = tempMat[2][2];

  R[6][6] = tempMat[0][0], R[6][7] = tempMat[0][1], R[6][8] = tempMat[0][2];
  R[7][6] = tempMat[1][0], R[7][7] = tempMat[1][1], R[7][8] = tempMat[1][2];
  R[8][6] = tempMat[2][0], R[8][7] = tempMat[2][1], R[8][8] = tempMat[2][2];

  return R;
}