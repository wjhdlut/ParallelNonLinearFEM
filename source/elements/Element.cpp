
#include <elements/Element.h>
#include <util/DataStructure.h>

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
    else{
      m_props[iter.key()] = iter.value();
    }
  }
}

Element::~Element()
{
}


void Element::AppendNodalOutput(const std::string&outputName, const Matrix&outMatrix)
{
  std::vector<double> outw(outMatrix.size(), 1.0);
  int numOfNode = GlobalData::GetInstance()->m_nodes->GetNumOfNodes();
  std::string outWeightName = outputName + "Weight";
  if(GlobalData::GetInstance()->m_outputData.count(outputName) == 0)
  {
    GlobalData::GetInstance()->m_outputData[outputName] = Math::MatrixZeros(numOfNode, outMatrix.size());
    GlobalData::GetInstance()->m_outputName.emplace_back(outputName);

    GlobalData::GetInstance()->m_outputData[outWeightName] = Math::MatrixZeros(numOfNode, 1);
    GlobalData::GetInstance()->m_outputName.emplace_back(outWeightName);
  }

  Matrix &outMatrix1 = GlobalData::GetInstance()->m_outputData[outputName];
  Matrix &outWeight = GlobalData::GetInstance()->m_outputData[outWeightName];

  if((outMatrix[0].size() != outMatrix1[0].size()) || outMatrix1.size() != m_nodes.size())
    throw "Appended output vector has incorrect size.";

  for(int row = 0; row < numOfNode; row++)
  {
    for(int line = 0; line < outMatrix[0].size(); line++)
    {
      outMatrix1[m_nodes[row]][line] += outMatrix[row][line];
    }
    outWeight[m_nodes[row]][0] += outw[row];
  }
}