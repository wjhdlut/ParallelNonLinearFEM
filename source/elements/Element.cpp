
#include <elements/Element.h>
#include <util/ShapeFunctions.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/DataStructure.h>
#include <util/ObjectFactory.h>

#include <iostream>

Element::Element(const std::vector<int> &elemNodes, const nlohmann::json &modelProps)
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
    throw "Appended output vector has incorrect size.";
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
  m_history = m_current;
  m_current.clear();

  if(m_mat != nullptr) m_mat->CommitHistory();
}

void Element::GetMassMatrix(std::shared_ptr<ElementData>&elemDat)
{
  std::string elemType = ShapeFunctions::GetElemType(elemDat->m_coords);
  if(m_reductedIntegration) order = -1;
  ShapeFunctions::GetIntegrationPoints(xi, weight, elemType, order, method);

  std::string elemName = elemType + "ShapeFunctions";
  std::shared_ptr<ElementShapeFunctions>res
               = ObjectFactory::CreateObject<ElementShapeFunctions>(elemName);
  if(nullptr == res) throw "Unkonwn type " + elemType;
  
  MatrixXd jac;
  double detJac;
  for(int i = 0; i < xi.rows(); i++)
  {
    res->GetShapeFunction(xi.row(i));
    
    // compute jacobian matrix
    jac = elemDat->m_coords.transpose() * res->pHpxi;
    detJac = jac.determinant();

    GetNMatrix(res->H);
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

void Element::ComputeElemTimeStep(const std::shared_ptr<ElementShapeFunctions> &res,
                                  const MatrixXd &elemNodeCoords,
                                  const VectorXd &elemNodeDisp,
                                  const double detJac)
{
  m_vol = res->ComputeElemTimeStep(m_dtK1, m_elemDistortion, elemNodeCoords,
                                   elemNodeDisp, detJac, m_waveSpeed);
}

void Element::HourGlassTech(std::shared_ptr<ElementData>&elemDat,
                            const VectorXd &elemNodeDisp,
                            const std::shared_ptr<ElementShapeFunctions> &res,
                            const MatrixXd &pHpX)
{
  if(m_props.contains("hourGlass")){
    const nlohmann::json &hourGlassPara = m_props.at("hourGlass");  
    double qh = hourGlassPara.at("para");
    double ah = qh * m_rho * pow(m_vol, 2./3.) / 4.;
    
    elemDat->m_fint -= ah * res->HourGlassTech(elemDat, elemNodeDisp, m_waveSpeed,
                                               hourGlassPara, pHpX);
  }
}

void Element::Initialize(const nlohmann::json &modelProps)
{
  for(auto iter = modelProps.begin(); iter != modelProps.end(); iter++)
  {
    if("material" == iter.key()){
      const nlohmann::json &matProps = iter.value();
      m_mat = std::make_shared<MaterialManager>(matProps);
      double E = m_mat->GetMaterialPara("E");
      double nu = m_mat->GetMaterialPara("nu");
      m_rho = m_mat->GetMaterialPara("rho");
      m_waveSpeed = sqrt(E*(1.-nu)/((1+nu)*(1-2.*nu)*m_rho));
    }
    else if("dim" == iter.key()){
      if("2D" == modelProps["dim"])
        m_dofType = {"u", "v"};
      if("3D" == modelProps["dim"])
        m_dofType = {"u", "v", "w"};
    }
    else if("reducedIntegration" == iter.key()){
      m_reductedIntegration = modelProps["reducedIntegration"];
    }
    else{
      m_props[iter.key()] = iter.value();
    }
  }
  
  if(!modelProps.contains("dim")) m_dofType = {"u", "v"};
}