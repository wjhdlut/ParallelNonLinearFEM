#include <materials/BaseMaterial.h>

BaseMaterial::BaseMaterial(const nlohmann::json &props) : m_props(props)
{}

BaseMaterial::~BaseMaterial()
{}

void BaseMaterial::SetIter(int iIter)
{
  m_iIter = iIter;
}

int BaseMaterial::SetHistoryParameter(const std::string &name, double value)
{
  if(-1 == m_iIter) 
  {
    m_initHistory[name] = value;
    return 0;
  }

  if(m_current.size() == m_iIter) m_current.emplace_back(m_initHistory);

  m_current[m_iIter].at(name) = value;

  return 0;
}

double BaseMaterial::GetHistoryParameter(const std::string&name)
{
  if(0 == m_history.size()) return m_initHistory[name];
  else return m_history[m_iIter].at(name);
}

void BaseMaterial::CommitHistory()
{
  m_history.clear();
  for(auto h : m_current)
    m_history.emplace_back(h);
}

std::vector<double> BaseMaterial::GetStress(const std::shared_ptr<Kinematics> &kin,
                                            const std::vector<double> &increDisp,
                                            const Matrix &dhpi){
  return Math::MatrixAMultVecB(m_D, kin->strain);
}

double BaseMaterial::SetMaterialParamter(const std::string &name)
{
  if(!m_props.contains(name)) return 0.;
  if(m_props.at(name).is_string()){
    std::string E = m_props.at(name);
    return std::stod(E);
  }
  else{
    return m_props.at(name);
  }
}
