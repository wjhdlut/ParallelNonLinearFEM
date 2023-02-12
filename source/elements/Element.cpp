
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
  if(std::find(GlobalData::GetInstance()->m_outputName.begin(), GlobalData::GetInstance()->m_outputName.end(),
    outputName) == GlobalData::GetInstance()->m_outputName.end())
  {
    GlobalData::GetInstance()->m_outputData[outputName] = Math::MatrixZeros(numOfNode, outMatrix[0].size());
    GlobalData::GetInstance()->m_outputName.emplace_back(outputName);

    GlobalData::GetInstance()->m_outputData[outWeightName] = Math::MatrixZeros(numOfNode, 1);
    GlobalData::GetInstance()->m_outputName.emplace_back(outWeightName);
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