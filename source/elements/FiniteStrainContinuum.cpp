
#include <elements/FiniteStrainContinuum.h>
#include <util/ShapeFunctions.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/Kinematics.h>
#include <iostream>

FiniteStrainContinuum::FiniteStrainContinuum(const std::vector<int> &elemNodes,
                                             const nlohmann::json &modelProps)
                      : Element(elemNodes, modelProps)
{
  // m_dofType = {"u", "v"};
  // m_rho = m_mat->GetMaterialRho();
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
  outputData.setZero(elemDat->m_coords.rows(), res->numOfStress);

  for(int iGaussPoint = 0; iGaussPoint < xi.rows(); iGaussPoint++){
    // compute shape functions
    res->GetShapeFunction(xi.row(iGaussPoint));

    // compute jacobian matrix
    jac = elemDat->m_coords.transpose() * res->pHpxi;
    // jac = Math::MatrixATransMultB(elemDat->m_coords, res->pHpxi);

    // Compute time step based on the deformation
    ComputeElemTimeStep(res, elemDat, jac.determinant());
    
    // compute the derivative of shape function about physical coordinate
    pHpX = res->pHpxi * jac.inverse();

    // compute deforamtion gradient
    GetKinematics(pHpX, elemDat->m_state);

    // compute tangent modulue matrix
    // D = m_mat->GetTangMatrix(kin->F);
    D = m_mat->GetTangMatrix();
    
    // compute strain matrix
    GetBMatrix(pHpX, kin->F);
    // kin->strain = B * elemDat->m_Dstate;
    
    // compute linear stiffness matrix
    elemDat->m_stiff += jac.determinant() * weight[iGaussPoint] * (B.transpose() * m_mat->GetTangMatrix() * B);

    // compute stress matrix
    sigma = m_mat->GetStress(kin, elemDat->m_Dstate);
    Stress2Matrix(sigma);
    
    // compute nonlinear strain matrix
    GetBNLMatrix(pHpX);

    // compute nonlinear stiffness matrix
    elemDat->m_stiff += jac.determinant() * weight[iGaussPoint] * (Bnl.transpose() * T * Bnl);

    // compute internal force vector
    elemDat->m_fint += jac.determinant() * weight[iGaussPoint] * (B.transpose() * sigma);
    
    // Hour-Glass method
    HourGlassTech(elemDat, res, pHpX);
    // std::cout << "xi = \n" << xi.row(iGaussPoint) << std::endl;
    // std::cout << "B = \n" << B << std::endl;
    // std::cout << "Bnl = \n" << Bnl << std::endl;
    // std::cout << "D = \n" << D << std::endl;
    // std::cout << "sigma = \n" << sigma << std::endl;
    // std::cout << "T = \n" << T<< std::endl;

    // compute output stress matrix
    outputData += Math::VecCross(VectorXd::Ones(elemDat->m_coords.rows()), sigma);
  }
  // std::cout << "elemDat->m_stiff\n" << elemDat->m_stiff << std::endl;
  // std::cout << "elemDat->m_fint\n" << elemDat->m_fint << std::endl;
  elemDat->m_outputData = 1./xi.rows() * outputData;
}

void FiniteStrainContinuum::GetKinematics(const MatrixXd &dphi,
                                          const VectorXd &elState)
{
  // compute deformation gradient tensor
  int numOfDim = dphi.cols();
  int numOfNode = dphi.rows();
  MatrixXd eleStateMat(numOfDim, numOfNode);
  eleStateMat = Math::ConvertVecToMat(numOfNode, numOfDim, elState);
  
  kin = std::make_shared<Kinematics>(numOfDim);
  kin->F += eleStateMat.transpose() * dphi;
  
  // compute Green-Lagrange strain tensor
  MatrixXd rightCauchyGreen = kin->F.transpose() * kin->F;
  kin->E = 0.5 * (rightCauchyGreen - MatrixXd::Identity(numOfDim, numOfDim));
  kin->SetStrainVector();
}


void FiniteStrainContinuum::GetBMatrix(const MatrixXd&dphi,
                                       const MatrixXd&F)
{
  int numOfDim = dphi.cols();
  int numOfNode = dphi.rows();

  if (2 == numOfDim)
  {
    B = MatrixXd::Zero(3, numOfDim * numOfNode);
    
    for (int i = 0; i < numOfNode; i++)
    {
      B(0, 2 * i + 0) = dphi(i, 0) * F(0, 0);
      B(0, 2 * i + 1) = dphi(i, 0) * F(1, 0);

      B(1, 2 * i + 0) = dphi(i, 1) * F(0, 1);
      B(1, 2 * i + 1) = dphi(i, 1) * F(1, 1);

      B(2, 2 * i + 0) = dphi(i, 1) * F(0, 0) + dphi(i, 0) * F(0, 1);
      B(2, 2 * i + 1) = dphi(i, 0) * F(1, 1) + dphi(i, 1) * F(1, 0);
    }
  }
  else if(3 == numOfDim)
  {
    B = MatrixXd::Zero(6, numOfDim * numOfNode);
    
    for(int i = 0; i < numOfNode; i++)
    {
      B(0, 3 * i + 0) = dphi(i, 0) * F(0, 0);
      B(0, 3 * i + 1) = dphi(i, 0) * F(1, 0);
      B(0, 3 * i + 2) = dphi(i, 0) * F(2, 0);

      B(1, 3 * i + 0) = dphi(i, 1) * F(0, 1);
      B(1, 3 * i + 1) = dphi(i, 1) * F(1, 1);
      B(1, 3 * i + 2) = dphi(i, 1) * F(2, 1);

      B(2, 3 * i + 0) = dphi(i, 2) * F(0, 2);
      B(2, 3 * i + 1) = dphi(i, 2) * F(1, 2);
      B(2, 3 * i + 2) = dphi(i, 2) * F(2, 2);

      B(3, 3 * i + 0) = dphi(i, 1) * F(0, 0) + dphi(i, 0) * F(0, 1);
      B(3, 3 * i + 1) = dphi(i, 0) * F(1, 1) + dphi(i, 1) * F(1, 0);
      B(3, 3 * i + 2) = dphi(i, 0) * F(2, 1) + dphi(i, 1) * F(2, 0);

      B(4, 3 * i + 0) = dphi(i, 1) * F(0, 2) + dphi(i, 2) * F(0, 1);
      B(4, 3 * i + 1) = dphi(i, 2) * F(1, 1) + dphi(i, 1) * F(1, 2);
      B(4, 3 * i + 2) = dphi(i, 2) * F(2, 1) + dphi(i, 1) * F(2, 2);

      B(5, 3 * i + 0) = dphi(i, 2) * F(0, 0) + dphi(i, 0) * F(0, 2);
      B(5, 3 * i + 1) = dphi(i, 0) * F(1, 2) + dphi(i, 2) * F(1, 0);
      B(5, 3 * i + 2) = dphi(i, 0) * F(2, 2) + dphi(i, 2) * F(2, 0);
    }
  }
}

void FiniteStrainContinuum::Stress2Matrix(const VectorXd&stress)
{
  if(3 == stress.size())
  {
    T = MatrixXd::Zero(4, 4);

    T(0, 0) = stress(0);
    T(1, 1) = stress(1);
    T(0, 1) = stress(2);
    T(1, 0) = stress(2);

    T(2, 2) = stress(0);
    T(3, 3) = stress(1);
    T(2, 3) = stress(2);
    T(3, 2) = stress(2);
  }
  else if(6 == stress.size())
  {
    T = MatrixXd::Zero(9, 9);

    T(0, 0) = stress(0), T(0, 1) = stress(3), T(0, 2) = stress(5);
    T(1, 0) = stress(3), T(1, 1) = stress(1), T(1, 2) = stress(4);
    T(2, 0) = stress(5), T(2, 1) = stress(4), T(2, 2) = stress(2);

    T(3, 3) = stress(0), T(3, 4) = stress(3), T(3, 5) = stress(5);
    T(4, 3) = stress(3), T(4, 4) = stress(1), T(4, 5) = stress(4);
    T(5, 3) = stress(5), T(5, 4) = stress(4), T(5, 5) = stress(2);

    T(6, 6) = stress(0), T(6, 7) = stress(3), T(6, 8) = stress(5);
    T(7, 6) = stress(3), T(7, 7) = stress(1), T(7, 8) = stress(4);
    T(8, 6) = stress(5), T(8, 7) = stress(4), T(8, 8) = stress(2);
  }
}

void FiniteStrainContinuum::GetBNLMatrix(const MatrixXd &dphi)
{
  if(dphi.rows() == 0) throw "the deratative of shape function to form BNL is empty";

  int numOfNode = dphi.rows();
  int numOfDim = dphi.cols();
  
  if(2 == numOfDim)
  {
    Bnl = MatrixXd::Zero(4, numOfNode*numOfDim);
    
    for (int i = 0; i < numOfNode; i++)
    {
      Bnl(0, numOfDim * i + 0) = dphi(i, 0);
      Bnl(1, numOfDim * i + 0) = dphi(i, 1);

      Bnl(2, numOfDim * i + 1) = dphi(i, 0);
      Bnl(3, numOfDim * i + 1) = dphi(i, 1);
    }
  }
  else if(3 == numOfDim)
  {
    Bnl = MatrixXd::Zero(9, numOfNode*numOfDim);
    
    for (int i = 0; i < numOfNode; i++)
    {
      Bnl(0, numOfDim * i + 0) = dphi(i, 0);
      Bnl(1, numOfDim * i + 0) = dphi(i, 1);
      Bnl(2, numOfDim * i + 0) = dphi(i, 2);

      Bnl(3, numOfDim * i + 1) = dphi(i, 0);
      Bnl(4, numOfDim * i + 1) = dphi(i, 1);
      Bnl(5, numOfDim * i + 1) = dphi(i, 2);

      Bnl(6, numOfDim * i + 2) = dphi(i, 0);
      Bnl(7, numOfDim * i + 2) = dphi(i, 1);
      Bnl(8, numOfDim * i + 2) = dphi(i, 2);
  }
  }
}

