#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>
#include <vector>

#include <nlohmann/json.hpp>
#include <materials/MaterialManager.h>

struct ElementData
{
  inline ElementData(std::vector<double>&elemState, std::vector<double>&elemDstate){
    m_state = elemState;
    m_Dstate = elemDstate;
    int nDof = elemState.size();
    for(int row = 0; row < nDof; row++){
      m_fint.emplace_back(0.);
      m_lumped.emplace_back(0.);
      std::vector<double> temp(nDof, 0.);
      m_stiff.emplace_back(temp);
      m_mass.emplace_back(temp);
    }
  }

  std::vector<double> m_state;
  std::vector<double> m_Dstate;
  std::vector<std::vector<double>> m_stiff;
  std::vector<double> m_fint;
  std::vector<std::vector<double>> m_mass;
  std::vector<double> m_lumped;
  std::vector<std::string> m_outLabel;
  std::vector<std::vector<double>> m_coords;
};

class Element
{
public:
  Element(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  virtual ~Element();

  inline std::vector<int> GetNodes(){
    return m_nodes;
  }

  inline std::vector<std::string> GetDofType(){
    return m_dofType;
  }

  inline void MatReset(){
    m_mat->Reset();
  }

  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) = 0;

protected:
  inline int DofCount(){
    return m_nodes.size() * m_dofType.size();
  }

protected:
  std::vector<std::string> m_dofType;
  std::shared_ptr<MaterialManager> m_mat;
  nlohmann::json m_props;
  std::vector<int> m_nodes;
};

#endif // ELEMENT_H
