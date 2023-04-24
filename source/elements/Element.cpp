
#include <elements/Element.h>
#include <util/ShapeFunctions.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/DataStructure.h>
#include <util/ObjectFactory.h>

Element::Element(const std::vector<int> &elemNodes, const nlohmann::json &modelProps)
        : m_nodes(elemNodes)
{
  for(auto iter = modelProps.begin(); iter != modelProps.end(); iter++)
  {
    if("material" == iter.key())
    {
      const nlohmann::json &matProps = iter.value();
      m_mat = std::make_shared<MaterialManager>(matProps);
    }
    else if("dim" == iter.key())
    {
      if("2D" == modelProps["dim"])
        m_dofType = {"u", "v"};
      if("3D" == modelProps["dim"])
        m_dofType = {"u", "v", "w"};
    }
    else if("reducedIntegration" == iter.key())
    {
      m_reductedIntegration = modelProps["reducedIntegration"];
    }
    else{
      m_props[iter.key()] = iter.value();
    }
  }
  if(!modelProps.contains("dim")) m_dofType = {"u", "v"};
  
  double E = m_mat->GetMaterialPara("E");
  double nu = m_mat->GetMaterialPara("nu");
  m_rho = m_mat->GetMaterialPara("rho");
  m_waveSpeed = sqrt(E*(1.-nu)/((1+nu)*(1-2.*nu)*m_rho));
}

Element::~Element()
{
}


void Element::AppendNodalOutput(const std::string&outputName, const Matrix&outMatrix)
{
  std::vector<double> outw(outMatrix.size(), 1.0);
  int numOfNode = GlobalData::GetInstance()->m_nodes->GetNumOfNodes();
  std::string outWeightName = outputName + "Weight";
  if(std::find(GlobalData::GetInstance()->m_outputName.begin(), GlobalData::GetInstance()->m_outputName.end(),
    outputName) == GlobalData::GetInstance()->m_outputName.end())
  {
    GlobalData::GetInstance()->m_outputData[outputName] = Math::MatrixZeros(numOfNode, outMatrix[0].size());
    GlobalData::GetInstance()->m_outputName.emplace_back(outputName);

    GlobalData::GetInstance()->m_outputData[outWeightName] = Math::MatrixZeros(numOfNode, 1);
    // GlobalData::GetInstance()->m_outputName.emplace_back(outWeightName);
  }

  Matrix &outMatrix1 = GlobalData::GetInstance()->m_outputData[outputName];
  Matrix &outWeight = GlobalData::GetInstance()->m_outputData[outWeightName];

  if((outMatrix[0].size() != outMatrix1[0].size()) || outMatrix.size() != m_nodes.size())
    throw "Appended output vector has incorrect size.";
  std::vector<int> index = GlobalData::GetInstance()->m_dofs->GetIndex(m_nodes);
  for(int row = 0; row < index.size(); row++)
  {
    for(int line = 0; line < outMatrix[0].size(); line++)
    {
      outMatrix1[index[row]][line] += outMatrix[row][line];
    }
    outWeight[index[row]][0] += outw[row];
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
  ShapeFunctions::GetIntegrationPoints(xi, weight, elemType, order, method);

  std::string elemName = elemType + "ShapeFunctions";
  std::shared_ptr<ElementShapeFunctions>res = ObjectFactory::CreateObject<ElementShapeFunctions>(elemName);
  if(nullptr == res) throw "Unkonwn type " + elemType;
  
  int count = 0;
  Matrix jac;
  double detJac;
  for(auto iXi : xi)
  {
    res->GetShapeFunction(iXi);
    
    // compute jacobian matrix
    jac = Math::MatrixATransMultB(elemDat->m_coords, res->pHpxi);
    detJac = Math::MatrixDet(jac);

    GetNMatrix(res->H);
    elemDat->m_mass = Math::MatrixAdd(m_rho*weight[count]*detJac, elemDat->m_mass,
                                      Math::MatrixATransMultB(N, N));
    count += 1;
  }
  
  // Compute lumped mass matrix
  for(int i = 0; i < elemDat->m_mass.size(); i++){
    elemDat->m_lumped[i] = Math::VecSum(elemDat->m_mass[i]);
  }
}

void Element::GetNMatrix(const std::vector<double> &H)
{
  int numOfDof = m_dofType.size();
  std::vector<double> temp(numOfDof*H.size());
  N.resize(numOfDof, temp);

  int lineIndex = 0;
  for(auto iH : H)
  {
    for(int i = 0; i < numOfDof; i++)
      N[i][lineIndex*numOfDof+i] = iH;
    
    lineIndex += 1;
  }
}