

#include <elements/SmallStrainContinuum.h>
#include <util/ShapeFunctions.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>

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

    // compute the derivative of shape function about physical coordinate
    invJac = Math::MatrixInverse(jac);
    pHpX = Math::MatrixAMultB(res->pHpxi, invJac);

    // compute strain matrix B
    B = GetBMatrix(pHpX);

    kin = GetKinematics(B, elemDat->m_state);

    // compute stress vector
    sigma = m_mat->GetStress(kin);
 
    // compute tangent modulue matrix
    D = m_mat->GetTangMatrix();

    // compute stiffness matrix
    detJac *= weight[count];
    tempMatrix = Math::MatrixATransMultB(B, Math::MatrixAMultB(m_mat->GetTangMatrix(), B));
    elemDat->m_stiff = Math::MatrixAdd(detJac, elemDat->m_stiff, tempMatrix);

    // compute internal force vector
    tempVec = Math::MatrixATransMultVecB(B, sigma);
    elemDat->m_fint = Math::VecAdd(detJac, elemDat->m_fint, tempVec);

    // compute output stress matrix
    tempMatrix = Math::VecOuter(std::vector<double>(elemDat->m_coords.size(), 1.0), sigma);
    outputData = Math::MatrixAdd(1., outputData, tempMatrix);

    count += 1;
  }
  elemDat->m_outputData = Math::MatrixScale(1./xi.size(), outputData);
}

std::shared_ptr<Kinematics> SmallStrainContinuum::GetKinematics(const Matrix &B, 
                                                                const std::vector<double> &elState)
{
  int numOfDim = 0;
  if(3 == B.size())
    numOfDim = 2;
  else if (6 == B.size())
    numOfDim = 3;
  else
    throw "Size of Strain Matrix B in SmallStrainContinuum::GetKinematics is Wrong";
  
  std::shared_ptr<Kinematics> kin = std::make_shared<Kinematics>(numOfDim);

  kin->strain = Math::MatrixAMultVecB(B, elState);

  return kin;
}

Matrix SmallStrainContinuum::GetBMatrix(const Matrix &dphi)
{
  int numOfDim = dphi[0].size();
  int numOfNode = dphi.size();
  std::vector<double> temp(numOfDim*numOfNode, 0.);
  Matrix B(3, temp);

  int count = 0;
  for(auto dp : dphi)
  {
    B[0][2*count + 0] = dp[0];
    B[1][2*count + 1] = dp[1];

    B[2][2*count + 0] = dp[1];
    B[2][2*count + 1] = dp[0];

    count += 1;
  }

  return B;
}