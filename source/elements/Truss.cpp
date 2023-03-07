#include <elements/Truss.h>
#include <util/Transformations.h>

Truss::Truss(const std::vector<int> &elemNode, const nlohmann::json &modelProps)
      : Element(elemNode, modelProps)
{
  SetHistoryParameter("sigma", 0.);
  
  CommitHistory();

  if(modelProps.at("E").is_string()){
    std::string E = modelProps.at("E");
    m_E = std::stod(E);
  }
  else{
    m_E = modelProps.at("E");
  }

  if(modelProps.at("Area").is_string()){
    std::string area = modelProps.at("Area");
    m_area = std::stod(area);
  }
  else{
    m_area = modelProps.at("Area");
  }
}

Truss::~Truss()
{}

void Truss::GetTangentStiffness(std::shared_ptr<ElementData> &elemDat)
{
  // Total displacement
  std::vector<double> a = Transformations::ToElementCoordinates(elemDat->m_state, elemDat->m_coords);
  
  // Increment displacement at current load step
  std::vector<double> Da = Transformations::ToElementCoordinates(elemDat->m_Dstate, elemDat->m_coords);
  
  // Converfenced displacement at last load step
  std::vector<double> a0 = Math::VecAdd(-1, a, Da);

  // Length of Truss element at the initial configuration
  m_l0 = Math::VecNorm(Math::VecAdd(-1., elemDat->m_coords[1], elemDat->m_coords[0]));

  double epsilon = 0., dEpsilon = 0.;
  GetStrain(epsilon, dEpsilon, a, a0);

  // Compute the stress increment (multiplied with the undeformed cross-sectional area)
  double dSigma = m_E * dEpsilon;

  // Compute the current stress (multiplied with the undeformed cross-sectional area)
  double sigma = GetHistoryParameter("sigma") + dSigma;

  // Update the history parameter
  SetHistoryParameter("sigma", sigma);

  // Compute BL in the element coordinate system
  std::vector<double> B = GetBMatrix(a);

  // Compute the element stiffness in the element coordinate system
  // KL = E*A0*l0*B*B^T
  Matrix KL = Math::MatrixScale(m_E * m_area * m_l0, Math::VecOuter(B, B));

  /**---------------------------------------------------------------------
   * KNL = A0*sigma/l0*[ 1,  0, -1,  0;
   *                     0,  1,, 0, -1;
   *                    -1,  0,  1,  0;
   *                     0, -1,  0,  1];
   * ---------------------------------------------------------------------*/
  Matrix KNL = GetNonLinearStiffMatrix(sigma, m_area);

  Matrix elemStiff = Math::MatrixAdd(1., KL, KNL);

  // Rotate element tangent stiffness to the global coordinate system
  elemDat->m_stiff = Transformations::ToGlobalCoordinates(elemStiff, elemDat->m_coords);
  
  // Compute the element internal force vector in the element coordinate system
  std::vector<double> elemFint = Math::VecScale(m_l0 * sigma * m_area, B);

  // Rotate element fint to the global coordinate system
  elemDat->m_fint = Transformations::ToGlobalCoordinates(elemFint, elemDat->m_coords);
}

void Truss::GetStrain(double &epsilon, double &dEpsilon,
                      const std::vector<double> &a, const std::vector<double> &a0)
{
  /* strain at the current load step
   * strain = (u2-u1)/l0 + 0.5*((u2-u1)/l0)^2 + 0.5*((v2-v1)/l0)^2;
   */
  epsilon = (a[2] - a[0])/m_l0 + 0.5 * pow((a[2]-a[0])/m_l0, 2) + 0.5 * pow((a[3]-a[1])/m_l0, 2);
  
  // Strain based on the displacement at the last convergenced load step
  double epsilon0 = (a0[2] - a0[0])/m_l0 + 0.5 * pow((a0[2]-a0[0])/m_l0, 2) + 0.5 * pow((a0[3]-a0[1])/m_l0, 2);

  // Compute the Strain Increment
  dEpsilon = epsilon - epsilon0;
}

std::vector<double> Truss::GetBMatrix(const std::vector<double> &a)
{
  /**
   * B = [-(1+(u2-u1)/l0), -(v2-v1)/l0, (1+(u2-u1)/l0),(v2-v1)/l0]/l0;
   * 
   */
  std::vector<double> B(4, 0.);
  B[0] = -1./m_l0 * (1.+(a[2]-a[0])/m_l0);
  B[1] = -1./m_l0 * (a[3]-a[1])/m_l0;
  B[2] = -B[0];
  B[3] = -B[1];

  return B;
}

std::vector<std::vector<double>> Truss::GetNonLinearStiffMatrix(const double &sigma,
                                                                const double &area)
{
  std::vector<double> tempVec(4, 0.);
  std::vector<std::vector<double>> KNL(4, tempVec);

  KNL[0][0] = sigma * area / m_l0;
  KNL[1][1] = sigma * area / m_l0;

  KNL[0][2] = -sigma * area / m_l0;
  KNL[1][3] = -sigma * area / m_l0;

  KNL[2][0] = KNL[0][2];
  KNL[3][1] = KNL[1][3];

  KNL[2][2] = KNL[0][0];
  KNL[3][3] = KNL[1][1];

  return KNL;
}