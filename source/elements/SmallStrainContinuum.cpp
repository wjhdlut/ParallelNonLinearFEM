

#include <elements/SmallStrainContinuum.h>
#include <util/ShapeFunctions.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>

#include <algorithm>

SmallStrainContinuum::SmallStrainContinuum(const std::vector<int> &elemNodes,
                                           const nlohmann::json &modelProps)
                                          : Element(elemNodes, modelProps)
{
}

SmallStrainContinuum::~SmallStrainContinuum()
{}

void SmallStrainContinuum::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{
  std::string elemType = ShapeFunctions::GetElemType(elemDat->m_coords);
  elemDat->m_outLabel.emplace_back("stresses");
  
  if(m_reductedIntegration) order = -1;
  ShapeFunctions::GetIntegrationPoints(xi, weight, elemType, order, method);

  std::string elemName = elemType + "ShapeFunctions";
  std::shared_ptr<ElementShapeFunctions> res = ObjectFactory::CreateObject<ElementShapeFunctions>(elemName);
  if(nullptr == res) throw "Unknown type " + elemType;
  outputData = Math::MatrixZeros(elemDat->m_coords.size(), res->numOfStress);

  int count = 0;
  Matrix tempMatrix;
  std::vector<double> tempVec;
  for(auto iXi : xi)
  {
    // compute shape functions
    res->GetShapeFunction(iXi);

    // compute jacobian matrix
    jac = Math::MatrixATransMultB(elemDat->m_coords, res->pHpxi);

    detJac = Math::MatrixDet(jac);

    ComputeElemTimeStep(res, elemDat);

    // compute the derivative of shape function about physical coordinate
    invJac = Math::MatrixInverse(jac);
    pHpX = Math::MatrixAMultB(res->pHpxi, invJac);

    // compute strain matrix B
    GetBMatrix(pHpX);

    GetKinematics(elemDat->m_state);

    // compute stress vector
    sigma = m_mat->GetStress(kin, elemDat->m_Dstate, pHpX);
 
    // compute tangent modulue matrix
    D = m_mat->GetTangMatrix();

    // compute stiffness matrix
    detJac *= weight[count];
    tempMatrix = Math::MatrixATransMultB(B, Math::MatrixAMultB(D, B));
    elemDat->m_stiff = Math::MatrixAdd(detJac, elemDat->m_stiff, tempMatrix);

    // compute internal force vector
    tempVec = Math::MatrixATransMultVecB(B, sigma);
    elemDat->m_fint = Math::VecAdd(detJac, elemDat->m_fint, tempVec);
    
    // Hour-Glass method
    HourGlassTech(elemDat, res);

    // compute output stress matrix
    tempMatrix = Math::VecOuter(std::vector<double>(elemDat->m_coords.size(), 1.0), sigma);
    outputData = Math::MatrixAdd(1., outputData, tempMatrix);

    count += 1;
  }
  elemDat->m_outputData = Math::MatrixScale(1./xi.size(), outputData);
}

void SmallStrainContinuum::GetKinematics(const std::vector<double> &elState)
{
  int numOfDim = 0;
  if(3 == B.size())
    numOfDim = 2;
  else if (6 == B.size())
    numOfDim = 3;
  else
    throw "Size of Strain Matrix B in SmallStrainContinuum::GetKinematics is Wrong";
  
  kin = std::make_shared<Kinematics>(numOfDim);

  kin->strain = Math::MatrixAMultVecB(B, elState);
}

void SmallStrainContinuum::GetBMatrix(const Matrix &dphi)
{
  int numOfDim = dphi[0].size();
  int numOfNode = dphi.size();
  std::vector<double> temp(numOfDim*numOfNode, 0.);
  if (2 == numOfDim)
  {
    B.resize(3, temp);

    int count = 0;
    for (auto dp : dphi)
    {
      B[0][2 * count + 0] = dp[0];
      B[1][2 * count + 1] = dp[1];

      B[2][2 * count + 0] = dp[1];
      B[2][2 * count + 1] = dp[0];

      count += 1;
    }
  }
  else if(3 == numOfDim){
    B.resize(6, temp);

    int count = 0;
    for(auto dp : dphi)
    {
      B[0][2 * count + 0] = dp[0];
      B[1][2 * count + 1] = dp[1];
      B[2][2 * count + 2] = dp[2];

      B[3][2 * count + 0] = dp[1], B[3][2 * count + 1] = dp[0];
      B[4][2 * count + 1] = dp[2], B[4][2 * count + 2] = dp[1];
      B[5][2 * count + 0] = dp[2], B[5][2 * count + 2] = dp[0];
    }
  }
}

void SmallStrainContinuum::ComputeElemTimeStep(const std::shared_ptr<ElementShapeFunctions> &res,
                                               const std::shared_ptr<ElementData> &elemDat)
{
  int k1, k2, k3, k4;
  double areal = 1.0e20, aream = 0.;
  double e, g, f, atest;
  double x13, x24, fs, ft;
  for(int iFace = 0; iFace < res->ReturnFaceMatrix()->size(); iFace++)
  {
      k1 = res->ReturnFaceMatrix()->at(iFace)[0] - 1;
      k2 = res->ReturnFaceMatrix()->at(iFace)[1] - 1;
      k3 = res->ReturnFaceMatrix()->at(iFace)[2] - 1;
      k4 = res->ReturnFaceMatrix()->at(iFace)[3] - 1;
    
    e = 0., f = 0., g = 0.;
    for(int iDof = 0; iDof < m_dofType.size(); iDof++){
      x13 = elemDat->m_coords[k3][iDof] - elemDat->m_coords[k1][iDof];
      x24 = elemDat->m_coords[k4][iDof] - elemDat->m_coords[k2][iDof];

      fs = x13 - x24;
      ft = x13 + x24;
      
      
      e += fs * fs;
      f += fs * ft;
      g += ft * ft;
    }
    atest = e*g - f*f;

    aream = std::max(atest, aream);
    areal = std::min(atest, areal);
  }
  m_vol = 8. * detJac;
  double at = areal / aream;
  double dt = 4 * m_vol / sqrt(aream);
  
  dt /= m_waveSpeed;
  
  m_dtK1 = std::min(m_dtK1, dt);
  m_elemDistortion = std::min(m_elemDistortion, at);
}

void SmallStrainContinuum::HourGlassTech(std::shared_ptr<ElementData> &elemDat,
                                         const std::shared_ptr<ElementShapeFunctions> &res)
{
  if(m_props.contains("hourGlass")){
    const nlohmann::json &hourGlassPara = m_props.at("hourGlass");
    std::vector<double> temp = res->HourGlassTech(elemDat, m_waveSpeed, hourGlassPara, pHpX);
  
    double qh = hourGlassPara.at("para");
    double ah = qh * m_rho * pow(m_vol, 2./3.) / 4.;

    temp = Math::VecScale(ah, temp);
    elemDat->m_fint = Math::VecAdd(1., elemDat->m_fint, temp);
  }
}