
#include <elements/FiniteStrainContinuum.h>
#include <util/ShapeFunctions.h>
#include <elements/ElementShapeFunctions.h>
#include <util/Math.h>
#include <util/Kinematics.h>

FiniteStrainContinuum::FiniteStrainContinuum(const std::vector<int> &elemNodes,
                                             const nlohmann::json &modelProps)
                      : Element(elemNodes, modelProps)
{
  m_dofType = {"u", "v"};
}

FiniteStrainContinuum::~FiniteStrainContinuum()
{}

void FiniteStrainContinuum::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{
  std::string elemType = ShapeFunctions::GetElemType(elemDat->m_coords);
  elemDat->m_outLabel.emplace_back("stresses");
  
  ShapeFunctions::GetIntegrationPoints(xi, weight, elemType, order, method);
  
  std::string elemName = elemType + "ShapeFunctions";
  std::shared_ptr<ElementShapeFunctions>res = ObjectFactory::CreateObject<ElementShapeFunctions>(elemName);
  if(nullptr == res) throw "Unkonwn type " + elemType;
  outputData = Math::MatrixZeros(elemDat->m_coords.size(), res->numOfStress);
  
  int count = 0;
  for(auto iXi : xi){
    res->GetShapeFunction(iXi);

    // compute jacobian matrix
    jac = Math::MatrixATransMultB(elemDat->m_coords, res->pHpxi);

    detJac = Math::MatrixDet(jac);
    
    // compute the derivative of shape function about physical coordinate
    invJac = Math::MatrixInverse(jac);
    pHpX = Math::MatrixAMultB(res->pHpxi, invJac);
    weighti = Math::MatrixDet(jac) * weight[count];

    // compute deforamtion gradient
    kin = GetKinematics(pHpX, elemDat->m_state);

    sigma = m_mat->GetStress(kin);
    
    // compute tangent modulue matrix
    D = m_mat->GetTangMatrix();
    
    // compute strain matrix
    B = GetBMatrix(pHpX, kin->F);
    
    // compute linear stiffness matrix
    detJac *= weight[count];
    tempMatrix = Math::MatrixATransMultB(B, Math::MatrixAMultB(m_mat->GetTangMatrix(), B));
    elemDat->m_stiff = Math::MatrixAdd(detJac, elemDat->m_stiff, tempMatrix);

    // compute stress matrix
    T = Stress2Matrix(sigma);
    
    // compute nonlinear strain matrix
    Bnl = GetBNLMatrix(pHpX);

    // compute nonlinear stiffness matrix
    tempMatrix = Math::MatrixATransMultB(Bnl, Math::MatrixAMultB(T, Bnl));
    elemDat->m_stiff = Math::MatrixAdd(detJac, elemDat->m_stiff, tempMatrix);

    // compute internal force vector
    tempVec = Math::MatrixTAMultVecB(B, sigma);
    elemDat->m_fint = Math::VecAdd(detJac, elemDat->m_fint, tempVec);
    
    // compute output stress matrix
    tempMatrix = Math::VecOuter(std::vector<double>(elemDat->m_coords.size(), 1.0), sigma);
    outputData = Math::MatrixAdd(1., outputData, tempMatrix);
    
    count += 1;
  }
elemDat->m_outputData = Math::MatrixScale(1./xi.size(), outputData);
}

std::shared_ptr<Kinematics> FiniteStrainContinuum::GetKinematics(const std::vector<std::vector<double>> &dphi,
                                                                 const std::vector<double>&elState)
{
  // compute deformation gradient tensor
  int numOfDim = dphi[0].size();
  int numOfNode = dphi.size();
  
  std::shared_ptr<Kinematics> kin = std::make_shared<Kinematics>(numOfDim);
  for(int i = 0; i < numOfNode; i++)
    for(int j = 0 ; j < numOfDim; j++)
      for(int k = 0; k < numOfDim; k++)
        kin->F[j][k] += dphi[i][k] * elState[numOfDim*i + j];
  
  // compute Green-Lagrange strain tensor
  std::vector<std::vector<double>> rightCauchyGreen = Math::MatrixATransMultB(kin->F, kin->F);
  kin->E = Math::MatrixAdd(-1., rightCauchyGreen, Math::MatrixEye(numOfDim));
  kin->E = Math::MatrixScale(0.5, kin->E);
  kin->SetStrainVector();

  return kin;
}


std::vector<std::vector<double>> FiniteStrainContinuum::GetBMatrix(const std::vector<std::vector<double>>&dphi,
                                                                   const std::vector<std::vector<double>>&F)
{
  int numOfDim = dphi[0].size();
  int numOfNode = dphi.size();
  std::vector<double> temp(numOfDim*numOfNode, 0.);
  std::vector<std::vector<double>> B(3, temp);
  int count = 0;
  for(auto dp : dphi)
  {
    B[0][2*count + 0] = dp[0] * F[0][0];
    B[0][2*count + 1] = dp[0] * F[1][0];

    B[1][2*count + 0] = dp[1] * F[0][1];
    B[1][2*count + 1] = dp[1] * F[1][1];

    B[2][2*count + 0] = dp[1] * F[0][0] + dp[0] * F[0][1];
    B[2][2*count + 1] = dp[0] * F[1][1] + dp[1] * F[1][0];

    count += 1;
  }

  return B;
}

std::vector<std::vector<double>> FiniteStrainContinuum::Stress2Matrix(const std::vector<double>&stress)
{
  std::vector<double> tempVec(4, 0.);
  std::vector<std::vector<double>> T(4, tempVec);

  T[0][0] = stress[0];
  T[1][1] = stress[1];
  T[0][1] = stress[2];
  T[1][0] = stress[2];

  T[2][2] = stress[0];
  T[3][3] = stress[1];
  T[2][3] = stress[2];
  T[3][2] = stress[2];
  return T;
}

Matrix FiniteStrainContinuum::GetBNLMatrix(const Matrix &dphi)
{
  if(dphi.size() == 0) throw "the deratative of shape function to form BNL is empty";

  int numOfNode = dphi.size();
  int numOfDim = dphi[0].size();
  std::vector<double> temp(numOfDim*numOfNode, 0.);
  Matrix Bnl(4, temp);

  int count = 0;
  for(auto dp : dphi){
    Bnl[0][numOfDim * count + 0] = dp[0];
    Bnl[1][numOfDim * count + 0] = dp[1];

    Bnl[2][numOfDim * count + 1] = dp[0];
    Bnl[3][numOfDim * count + 1] = dp[1];

    count += 1;
  }

  return Bnl;
}