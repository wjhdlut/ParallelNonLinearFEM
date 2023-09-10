

#include <elements/SmallStrainContinuum.h>
#include <util/ShapeFunctions.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>

#include <algorithm>
#include <iostream>

SmallStrainContinuum::SmallStrainContinuum(const std::vector<int> &elemNodes,
                                           const nlohmann::json &modelProps)
                                          : Element(elemNodes, modelProps)
{
}

SmallStrainContinuum::~SmallStrainContinuum()
{
}

void SmallStrainContinuum::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{
  std::string elemType = ShapeFunctions::GetElemType(elemDat->m_coords);
  elemDat->m_outLabel.emplace_back("stresses");
  
  if(m_reductedIntegration) order = -1;
  ShapeFunctions::GetIntegrationPoints(xi, weight, elemType, order, method);

  std::string elemName = elemType + "ShapeFunctions";
  std::shared_ptr<ElementShapeFunctions> res = ObjectFactory::CreateObject<ElementShapeFunctions>(elemName);
  if(nullptr == res) throw "Unknown type " + elemType;
  outputData = MatrixXd::Zero(elemDat->m_coords.rows(), res->numOfStress);
  
  for(int i = 0; i < xi.rows(); i++)
  {
    // compute shape functions
    res->GetShapeFunction(xi.row(i));

    // compute jacobian matrix
    jac = elemDat->m_coords.transpose() * res->pHpxi;

    // Compute time step based on the deformation
    ComputeElemTimeStep(res, elemDat->m_coords, 
                        VectorXd::Zero(elemDat->m_state.size()), jac.determinant());

    // compute the derivative of shape function about physical coordinate
    pHpX = res->pHpxi * jac.inverse();

    // compute strain matrix B
    GetBMatrix(pHpX);

    GetKinematics(elemDat->m_state);

    // compute stress vector
    sigma = m_mat->GetStress(kin, elemDat->m_Dstate, pHpX);
 
    // compute tangent modulue matrix
    D = m_mat->GetTangMatrix();

    // compute stiffness matrix
    elemDat->m_stiff += jac.determinant() * weight(i) * B.transpose() * D * B;

    // compute internal force vector
    elemDat->m_fint += jac.determinant() * weight(i) * B.transpose() * sigma;
    
    // Hour-Glass method
    HourGlassTech(elemDat, VectorXd::Zero(elemDat->m_state.size()), res, pHpX);

    // compute output stress matrix
    outputData += Math::VecCross(VectorXd::Ones(elemDat->m_coords.rows()), sigma);
  }
  elemDat->m_outputData = 1./xi.rows() * outputData;
}

void SmallStrainContinuum::GetKinematics(const VectorXd &elState)
{
  int numOfDim = 0;
  if(3 == B.rows())
    numOfDim = 2;
  else if (6 == B.rows())
    numOfDim = 3;
  else
    throw "Size of Strain Matrix B in SmallStrainContinuum::GetKinematics is Wrong";
  
  kin = std::make_shared<Kinematics>(numOfDim);
  kin->strain = B * elState;
}

void SmallStrainContinuum::GetBMatrix(const MatrixXd &dphi)
{
  int numOfDim = dphi.cols();
  int numOfNode = dphi.rows();
  
  if (2 == numOfDim)
  {
    B = MatrixXd::Zero(3, numOfDim*numOfNode);
    
    for(int i = 0; i < numOfNode; i++)
    {
      B(0, numOfDim * i + 0) = dphi(i, 0);
      B(1, numOfDim * i + 1) = dphi(i, 1);

      B(2, numOfDim * i + 0) = dphi(i, 1);
      B(2, numOfDim * i + 1) = dphi(i, 0);
    }
  }
  else if(3 == numOfDim){
    B = MatrixXd::Zero(6, numOfDim*numOfNode);

    for(int i = 0; i < numOfNode; i++)
    {
      B(0, numOfDim * i + 0) = dphi(i, 0);
      B(1, numOfDim * i + 1) = dphi(i, 1);
      B(2, numOfDim * i + 2) = dphi(i, 2);

      B(3, numOfDim * i + 0) = dphi(i, 1), B(3, numOfDim * i + 1) = dphi(i, 0);
      B(4, numOfDim * i + 1) = dphi(i, 2), B(4, numOfDim * i + 2) = dphi(i, 1);
      B(5, numOfDim * i + 0) = dphi(i, 2), B(5, numOfDim * i + 2) = dphi(i, 0);
    }
  }
}