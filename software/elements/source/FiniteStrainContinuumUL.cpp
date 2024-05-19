/**
 * @File Name:     FiniteStrainContinuumUL.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/FiniteStrainContinuumUL.h"
#include "../shapefunctions/include/ElementShapeFunctions.h"
#include "../../util/include/ShapeFunctions.h"
#include "../../util/include/Kinematics.h"

#include <iostream>

FiniteStrainContinuumUL::FiniteStrainContinuumUL(const std::string &elemShape,
                                                 const std::vector<int> &elemNodes,
                                                 const nlohmann::json &modelProps)
                       : Element(elemNodes, modelProps)
{
  SetDofType(elemShape);

  CompGaussPointCoord(elemShape);

  InitializeHistoryVariables();
}

FiniteStrainContinuumUL::~FiniteStrainContinuumUL()
{}

void FiniteStrainContinuumUL::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{
  InsertElemOutputData(elemDat->m_outputData, "stresses");
  
  outputData.setZero(elemDat->m_coords.rows(), m_elemShapePtr->numOfStress);

  for(int iGaussPoint = 0; iGaussPoint < xi.rows(); iGaussPoint++){
    // compute shape functions
    m_elemShapePtr->GetShapeFunction(xi.row(iGaussPoint));

    nodeCoord = UpdateNodeCoords(elemDat);

    // compute jacobian matrix
    Jac = elemDat->m_coords.transpose() * m_elemShapePtr->pHpxi;
    jac = nodeCoord.transpose() * m_elemShapePtr->pHpxi;

    // Compute time step based on the deformation
    ComputeElemTimeStep(elemDat->m_coords, elemDat->m_state, jac.determinant());
    pHpX = m_elemShapePtr->pHpxi * Jac.inverse();

    // compute the derivative of shape function about 
    // physical coordinate in current configuration
    pHpx = m_elemShapePtr->pHpxi * jac.inverse();

    // compute deforamtion gradient
    GetKinematics(pHpX, elemDat);

    // compute tangent modulue matrix
    D = m_mat->GetTangMatrix();
    
    // compute strain matrix
    GetBMatrix(pHpx);
    if("AxiSymmetry" == m_analyseType)
      GetBMatrixForAxiSymmetry(B, m_elemShapePtr->H, nodeCoord);
    kin->strain = B * elemDat->m_state;
    
    // compute linear stiffness matrix
    elemDat->m_stiff += jac.determinant() * weight[iGaussPoint]
                      * (B.transpose() * m_mat->GetTangMatrix() * B);

    // compute stress matrix
    sigma = m_mat->GetStress(kin);
    Stress2Matrix(sigma);
    
    // compute nonlinear strain matrix
    GetBNLMatrix(pHpx);

    // compute nonlinear stiffness matrix
    elemDat->m_stiff += jac.determinant() * weight[iGaussPoint] * (Bnl.transpose() * T * Bnl);

    // compute internal force vector
    elemDat->m_fint += jac.determinant() * weight[iGaussPoint] * (B.transpose() * sigma);
    
    // Hour-Glass method
    HourGlassTech(elemDat, elemDat->m_state, pHpx);

    // compute output stress matrix
    outputData += Math::VecCross(VectorXd::Ones(elemDat->m_coords.rows()), sigma);
  }
  InsertElemOutputData(elemDat->m_outputData, "stresses", 1./xi.rows() * outputData);
}

void FiniteStrainContinuumUL::GetKinematics(const MatrixXd &dphi,
                                            const std::shared_ptr<ElementData> &elemDat)
{
  // compute deformation gradient tensor
  int numOfDim = dphi.cols();
  int numOfNode = dphi.rows();
  MatrixXd eleStateMat = Map<const MatrixXd>(elemDat->m_state.data(), numOfDim, numOfNode);
  
  kin = std::make_shared<Kinematics>(numOfDim);
  kin->F += eleStateMat * dphi;
  
  // compute Green-Lagrange strain tensor
  MatrixXd rightCauchyGreen = kin->F.transpose() * kin->F;
  kin->E = 0.5 * (rightCauchyGreen - MatrixXd::Identity(numOfDim, numOfDim));
  kin->SetStrainVector(kin->E);
}


void FiniteStrainContinuumUL::GetBMatrix(const MatrixXd&dphi)
{
  int numOfDim = dphi.cols();
  int numOfNode = dphi.rows();

  if (2 == numOfDim)
  {
    B = MatrixXd::Zero(3, numOfDim * numOfNode);
    
    for (int i = 0; i < numOfNode; i++)
    {
      B(0, 2 * i + 0) = dphi(i, 0);
      
      B(1, 2 * i + 1) = dphi(i, 1);

      B(2, 2 * i + 0) = dphi(i, 1);
      B(2, 2 * i + 1) = dphi(i, 0);
    }
  }
  else if(3 == numOfDim)
  {
    B = MatrixXd::Zero(6, numOfDim * numOfNode);
    
    for(int i = 0; i < numOfNode; i++)
    {
      B(0, 3 * i + 0) = dphi(i, 0);
      B(1, 3 * i + 1) = dphi(i, 1);
      B(2, 3 * i + 2) = dphi(i, 2);

      B(3, 3 * i + 0) = dphi(i, 1);
      B(3, 3 * i + 1) = dphi(i, 0);
      
      B(4, 3 * i + 1) = dphi(i, 2);
      B(4, 3 * i + 2) = dphi(i, 1);

      B(5, 3 * i + 0) = dphi(i, 2);
      B(5, 3 * i + 2) = dphi(i, 0);
    }
  }
}

void FiniteStrainContinuumUL::Stress2Matrix(const VectorXd&stress)
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

void FiniteStrainContinuumUL::GetBNLMatrix(const MatrixXd &dphi)
{
  if(dphi.rows() == 0){
    std::cout << "Catch Exception: "
              << "the deratative of shape function to form BNL is empty"
              << std::endl;
    exit(-1);
  }

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

MatrixXd FiniteStrainContinuumUL::UpdateNodeCoords(std::shared_ptr<ElementData>&elemDat)
{
  int numOfNode = elemDat->m_coords.rows();
  int numOfDim  = elemDat->m_coords.cols();
  MatrixXd nodeDisp = Map<MatrixXd>(elemDat->m_state.data(), numOfDim, numOfNode);
  
  return elemDat->m_coords + nodeDisp.transpose();
}