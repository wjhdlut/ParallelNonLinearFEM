/**
 * @File Name:     Element.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/Element.h"
#include "../../util/include/ShapeFunctions.h"
#include "../shapefunctions/include/ElementShapeFunctions.h"
#include "../../util/include/DataStructure.h"
#include "../../util/include/ObjectFactory.h"

#include <iostream>

Element::Element(const std::vector<int> &elemNodes,
                 const nlohmann::json &modelProps)
        : m_nodes(elemNodes)
{
  Initialize(modelProps);
}

Element::~Element()
{
}

void Element::AppendNodalOutput(const std::string&outputName, const MatrixXd&outMatrix)
{
  std::vector<double> outw(outMatrix.size(), 1.0);
  int numOfNode = GlobalData::GetInstance()->m_nodes->GetNumOfNodes();
  std::string outWeightName = outputName + "Weight";
  if(std::find(GlobalData::GetInstance()->m_outputName.begin(), GlobalData::GetInstance()->m_outputName.end(),
    outputName) == GlobalData::GetInstance()->m_outputName.end())
  {
    GlobalData::GetInstance()->m_outputData[outputName] = MatrixXd::Zero(numOfNode, outMatrix.cols());
    
    GlobalData::GetInstance()->m_outputName.emplace_back(outputName);

    GlobalData::GetInstance()->m_outputData[outWeightName] = MatrixXd::Zero(numOfNode, 1);
    // GlobalData::GetInstance()->m_outputName.emplace_back(outWeightName);
  }

  MatrixXd &outMatrix1 = GlobalData::GetInstance()->m_outputData[outputName];
  MatrixXd &outWeight = GlobalData::GetInstance()->m_outputData[outWeightName];

  if((outMatrix.cols() != outMatrix1.cols()) || outMatrix.rows() != m_nodes.size())
  {
    std::cout << "Catch Exception: "
              << "Appended output vector has incorrect size."
              << std::endl;
    exit(-1);
  }
  std::vector<int> index = GlobalData::GetInstance()->m_dofs->GetIndex(m_nodes);
  for(int row = 0; row < index.size(); row++)
  {
    for(int line = 0; line < outMatrix.cols(); line++)
    {
      outMatrix1(index[row], line) += outMatrix(row, line);
    }
    outWeight(index[row], 0) += outw[row];
  }
}

void Element::CommitHistory()
{
  if(0 != m_current.size()) m_history = m_current;
  m_current.clear();

  if(m_mat != nullptr) m_mat->CommitHistory();
}

void Element::GetMassMatrix(std::shared_ptr<ElementData>&elemDat)
{
  MatrixXd jac;
  double detJac;
  for(int i = 0; i < xi.rows(); i++)
  {
    m_elemShapePtr->GetShapeFunction(xi.row(i));
    
    // compute jacobian matrix
    jac = elemDat->m_coords.transpose() * m_elemShapePtr->pHpxi;
    detJac = jac.determinant();

    GetNMatrix(m_elemShapePtr->H);
    elemDat->m_mass += m_rho*weight[i]*detJac * (N.transpose() * N);
  }
  
  // Compute lumped mass matrix
  for(int i = 0; i < elemDat->m_mass.rows(); i++){
    elemDat->m_lumped[i] = elemDat->m_mass.row(i).sum();
  }
}

void Element::GetNMatrix(const VectorXd &H)
{
  int numOfDof = m_dofType.size();
  N = MatrixXd::Zero(numOfDof, numOfDof*H.size());

  for(int lineIndex = 0; lineIndex < H.size(); lineIndex++)
  {
    for(int i = 0; i < numOfDof; i++)
      N(i, lineIndex*numOfDof+i) = H(i);
  }
}

void Element::ComputeElemTimeStep(const MatrixXd &elemNodeCoords,
                                  const VectorXd &elemNodeDisp,
                                  const double detJac)
{
  m_vol = m_elemShapePtr->ComputeElemTimeStep(m_dtK1, m_elemDistortion, elemNodeCoords,
                                              elemNodeDisp, detJac, m_waveSpeed);
}

void Element::HourGlassTech(std::shared_ptr<ElementData>&elemDat,
                            const VectorXd &elemNodeDisp,
                            const MatrixXd &pHpX)
{
  if(m_props.contains("hourGlass")){
    const nlohmann::json &hourGlassPara = m_props.at("hourGlass");  
    double qh = hourGlassPara.at("para");
    double ah = qh * m_rho * pow(m_vol, 2./3.) / 4.;
    
    elemDat->m_fint -= ah * m_elemShapePtr->HourGlassTech(elemDat, elemNodeDisp, m_waveSpeed,
                                                          hourGlassPara, pHpX);
  }
}

void Element::Initialize(const nlohmann::json &modelProps)
{
  for(auto iter = modelProps.begin(); iter != modelProps.end(); iter++)
  {
    if("material" == iter.key()){
      SetMaterialParameters(modelProps);
    }
    else if("reducedIntegration" == iter.key()){
      m_reductedIntegration = modelProps["reducedIntegration"];
    }
    else if("analyseType" == iter.key()){
      SetAnalyseType(modelProps);
    }
    else{
      m_props[iter.key()] = iter.value();
    }
  }
}

void Element::SetMaterialParameters(const nlohmann::json &modelProps)
{
  nlohmann::json matProps = modelProps["material"];
  if (modelProps.contains("analyseType"))
    matProps["analyseType"] = modelProps["analyseType"];
  m_mat = std::make_shared<MaterialManager>(matProps);
  double E = m_mat->GetMaterialPara("E");
  double nu = m_mat->GetMaterialPara("nu");
  m_rho = m_mat->GetMaterialPara("rho");
  m_waveSpeed = sqrt(E * (1. - nu) / ((1 + nu) * (1 - 2. * nu) * m_rho));
}

void Element::SetAnalyseType(const nlohmann::json &modelProps)
{
  // For Axis Aymmetric Problem
  if ("AxiSymmetry" == modelProps["analyseType"])
  {
    m_analyseType = "AxiSymmetry";
    if (!modelProps.contains("axiSymmetry"))
    {
      std::cout << "Catch Exception: "
                << "Please Assign Axis of Symmetry!!!"
                << std::endl;
      exit(-1);
    }
    m_axiSymmetry = modelProps["axiSymmetry"];
  }
  // For Plane Stress Problem
  else if ("PlaneStrain" == modelProps["analyseType"])
  {
    m_analyseType = "PlaneStress";
  }
  // For Plane Strain Problem
  else if ("PlaneStrain" == modelProps["analyseType"])
  {
    m_analyseType = "PlaneStrain";
  }
}

void Element::SetDofType(const std::string &elemShape)
{
  std::string elemType = elemShape + "ShapeFunctions";
  m_elemShapePtr = ObjectFactory::CreateObject<ElementShapeFunctions>(elemType);
  if(nullptr == m_elemShapePtr){
    std::cout << "Catch Exception: "
              << "Unkonwn type " + elemType
              << std::endl;
    exit(-1);
  }
  m_dofType = m_elemShapePtr->ReturnDofType();
}

void Element::CompGaussPointCoord(const std::string &elemShape)
{
  if(m_reductedIntegration) order = -1;
  ShapeFunctions::GetIntegrationPoints(xi, weight, elemShape, order);
}

void Element::InitializeHistoryVariables()
{
  InitializeStress();
  
  InitializeStateVariable();

  CommitHistory();
}

void Element::InitializeStress()
{
  VectorXd tempVec = VectorXd::Zero(m_elemShapePtr->numOfStress * weight.size());
  SetHistoryParameter("sigma", tempVec);
}

void Element::InitializeStateVariable()
{
  VectorXd tempVec = VectorXd::Zero(m_elemShapePtr->numOfStress);
  SetHistoryParameter("stateVariable", tempVec);
}

std::vector<int> Element::CheckNodeBoundary(std::vector<int> &nodeChk,
                                            const std::unordered_map<int, std::vector<double>> &nodeForcePres)
{
  const std::unordered_map<int, std::vector<int>> &elemNodeOrdered = m_elemShapePtr->ReturnElemNodeOrdered();
  
  for(auto iter : nodeForcePres){
    int index = std::distance(m_nodes.begin(),
                std::find(m_nodes.begin(), m_nodes.end(), iter.first));
    nodeChk[index] = 1;
  }

  bool found = false;
  for(int iEdge = 1; iEdge <= elemNodeOrdered.size(); iEdge++)
  {
    for(int iNode = 0; iNode < m_nodes.size(); iNode++)
    {
      if((0 != nodeChk[iNode] && 0 == elemNodeOrdered.at(iEdge)[iNode])
      || (0 == nodeChk[iNode] && 0 != elemNodeOrdered.at(iEdge)[iNode])){
        continue;
      }
    }
    // Found the index of edge
    for(int iNode = 0; iNode < m_nodes.size(); iNode++){
      if(0 != elemNodeOrdered.at(iEdge)[iNode])
        nodeChk[elemNodeOrdered.at(iEdge)[iNode] - 1] = iNode + 1;
    }
  }
  std::vector<int> nodeAux;
  for(auto iter : nodeForcePres){
    for(int index = 0; index < nodeForcePres.size(); index++)
    {
      if (iter.first == m_nodes[nodeChk[index] - 1])
        nodeAux.emplace_back(iter.first);
    }
  }
  return nodeAux;
}

VectorXd Element::CompEquivalentNodeForce(const std::unordered_map<int, VectorXd> &nodeCoord,
                                          const std::unordered_map<int, std::vector<double>> &nodeForcePres)
{
  std::vector<int> nodeChk(m_nodes.size(), 0);
  VectorXd equivalentForce = VectorXd::Zero(DofCount());
  std::vector<int> nodeAux = CheckNodeBoundary(nodeChk, nodeForcePres);
  
  // Get Boundary Integration Point Data
  VectorXd boundaryWeight, boundaryH;
  MatrixXd boundaryXi, pBoundaryHpXi;
  m_elemShapePtr->GetBoundaryIntegrationPoint(boundaryXi, boundaryWeight);


  for(int iGaussPoint = 0; iGaussPoint < boundaryWeight.size(); iGaussPoint++){
    m_elemShapePtr->GetBoundaryShapeFunction(boundaryH, pBoundaryHpXi, boundaryXi.row(iGaussPoint));

    /*********************************************************************************
     * 
     * PGASH = [H] * Press stands for the values of Edge Load at gauss point iGaussPoint
     * DGASH = [pX/pxi pY/pxi] = [cons(thea), sin(thea)] stands for the derection
     * 
     *********************************************************************************/
    std::vector<double> PGASH(nodeForcePres.begin()->second.size(), 0.);
    std::vector<double> DGASH(PGASH);
    double cosThea = 0., sinThea = 0.;
    for(auto iter : nodeForcePres)
    {
      int ii = std::distance(nodeAux.begin(), std::find(nodeAux.begin(), nodeAux.end(), iter.first));
      for(int iDof = 0; iDof < iter.second.size(); iDof++)
      {
        PGASH[iDof] += iter.second[iDof] * boundaryH(ii);
        DGASH[iDof] += nodeCoord.at(iter.first)(iDof) * pBoundaryHpXi(ii, 0);
      }
    }
    double pX = DGASH[0] * PGASH[1] - DGASH[1] * PGASH[0];
    double pY = DGASH[0] * PGASH[0] + DGASH[1] * PGASH[1];

    for(int i = 0; i < nodeForcePres.size(); i++)
    {
      int iNode = nodeChk[i];
      int index = (iNode - 1) * m_dofType.size();
      equivalentForce(index + 0) += boundaryH(i) * pX * boundaryWeight(iGaussPoint);
      equivalentForce(index + 1) += boundaryH(i) * pY * boundaryWeight(iGaussPoint);
    }
  }
  return equivalentForce;
}

void Element::InsertElemOutputData(std::unordered_map<std::string, MatrixXd> &elemOutData,
                                   const std::string &name,
                                   const MatrixXd &outputData)
{
  if(0 == elemOutData.count(name))
    elemOutData.insert(std::pair<std::string, MatrixXd>(name, outputData));
  else
    elemOutData[name] = outputData;
}