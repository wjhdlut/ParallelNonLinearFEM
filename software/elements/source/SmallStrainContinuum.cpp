/**
 * @File Name:     SmallStrainContinuum.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <algorithm>
#include <iostream>

#include "../include/SmallStrainContinuum.h"
#include "../shapefunctions/include/ElementShapeFunctions.h"
#include "../../util/include/ShapeFunctions.h"

SmallStrainContinuum::SmallStrainContinuum(const std::string &elemShape,
                                           const std::vector<int> &elemNodes,
                                           const nlohmann::json &modelProps)
                     : Element(elemNodes, modelProps)
{
  Initialize(elemShape);
}

SmallStrainContinuum::~SmallStrainContinuum()
{
}

void SmallStrainContinuum::Initialize(const std::string &elemShape)
{
  SetDofType(elemShape);

  CompGaussPointCoord(elemShape);

  InitializeHistoryVariables();
}

void SmallStrainContinuum::GetTangentStiffness(std::shared_ptr<ElementData>&elemDat)
{
  InsertElemOutputData(elemDat->m_outputData, "stresses");
  
  if("AxiSymmetry" == m_analyseType)
    outputData = MatrixXd::Zero(elemDat->m_coords.rows(), 4);
  else
    outputData = MatrixXd::Zero(elemDat->m_coords.rows(), m_elemShapePtr->numOfStress);
  
  double tempWeight = 1.;
  // std::cout << "elemNodeCoord = \n" << elemDat->m_coords << std::endl;
  for(int i = 0; i < xi.rows(); i++)
  {
    // std::cout << i << "-th Gauss Point, " << xi.row(i) << std::endl;
    // compute shape functions
    m_elemShapePtr->GetShapeFunction(xi.row(i));

    // compute jacobian matrix
    jac = elemDat->m_coords.transpose() * m_elemShapePtr->pHpxi;

    // Compute time step based on the deformation
    ComputeElemTimeStep(elemDat->m_coords, 
                        VectorXd::Zero(elemDat->m_state.size()), jac.determinant());

    // compute the derivative of shape function about physical coordinate
    pHpX = m_elemShapePtr->pHpxi * jac.inverse();

    // compute strain matrix B
    GetBMatrix(pHpX);
    if("AxiSymmetry" == m_analyseType)
      tempWeight = GetBMatrixForAxiSymmetry(B, m_elemShapePtr->H, elemDat->m_coords);
    // std::cout << "B = \n" << B << std::endl;

    GetKinematics(elemDat);

    // compute stress vector
    sigma = m_mat->GetStress(kin);
    // std::cout << "sigma = " << sigma.transpose() << std::endl;
 
    // compute tangent modulue matrix
    D = m_mat->GetTangMatrix();
    Eigen::IOFormat fmt(10);
    // std::cout << "m_D = \n" << D.format(fmt) << std::endl;

    // compute stiffness matrix
    elemDat->m_stiff += jac.determinant() * weight(i) * B.transpose() * D * B * tempWeight;

    // compute internal force vector
    elemDat->m_fint += jac.determinant() * weight(i) * B.transpose() * sigma * tempWeight;
    
    // Hour-Glass method
    HourGlassTech(elemDat, VectorXd::Zero(elemDat->m_state.size()), pHpX);

    // compute output stress matrix
    outputData += Math::VecCross(VectorXd::Ones(elemDat->m_coords.rows()), sigma);
  }
  Eigen::IOFormat fmt(10);
  // std::cout << "elemStiff = \n" << elemDat->m_stiff.format(fmt) << std::endl;
  InsertElemOutputData(elemDat->m_outputData, "stresses", 1./xi.rows() * outputData);
}

void SmallStrainContinuum::GetKinematics(const std::shared_ptr<ElementData> &elemDat)
{
  int numOfDim = 0;
  if("PlaneStrain" == m_analyseType)
    numOfDim = 2;
  else if("PlaneStress" == m_analyseType)
    numOfDim = 2;
  else if("AxiSymmetry" == m_analyseType)
    numOfDim = 2;
  else if (6 == B.rows())
    numOfDim = 3;
  else{
    std::cout << "Catch Exception: "
              << "Size of Strain Matrix B in SmallStrainContinuum::GetKinematics is Wrong"
              << std::endl;
    exit(-1);
  }
  
  kin = std::make_shared<Kinematics>(numOfDim);
  if("AxiSymmetry" == m_analyseType){
    kin->strain.resize(4);
    kin->incremStrain.resize(4);
  }
  kin->strain = B * elemDat->m_state;
  kin->incremStrain = B * elemDat->m_Dstate;
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