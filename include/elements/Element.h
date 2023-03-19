#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>
#include <vector>

#include <nlohmann/json.hpp>
#include <materials/MaterialManager.h>

/**
 * @Brief:  element data data structure
 * 
 */
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

  std::vector<double> m_state;                   // element state variable vector such as displacement
  std::vector<double> m_Dstate;                  // increment of element state variable vector
  std::vector<std::vector<double>> m_stiff;      // element stiffness matrix
  std::vector<double> m_fint;                    // element internal force vector
  std::vector<std::vector<double>> m_mass;       // element mass matrix
  std::vector<double> m_lumped;                  // 
  std::vector<std::string> m_outLabel;           // output variable name
  std::vector<std::vector<double>> m_coords;     // element node coordinates
  std::vector<std::vector<double>> m_outputData; // output variable data
};

/**
 * @Brief: father class of element
 * 
 */
class Element
{
public:
  Element(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  
  virtual ~Element();

  /**
   * @Brief: Get element node index vector
   * 
   * @return std::vector<int>           element node index vector
   */
  inline std::vector<int> GetNodes(){
    return m_nodes;
  }
  
  /**
   * @Brief:  Get the Dof Type object
   * 
   * @return std::vector<std::string>   element dof type
   */
  inline std::vector<std::string> GetDofType(){
    return m_dofType;
  }

  
  inline void MatReset(){
    if(m_mat != nullptr)
      m_mat->Reset();
  }

  inline std::vector<double> GetHistoryParameter(const std::string&name){
    return m_history[name];
  }
  
  inline void SetHistoryParameter(const std::string&name, const std::vector<double> &value){
    m_current[name] = value;
  }

  /**
   * @Brief: Compute the Element Tangent Stiffness Matrix
   * 
   * @param elemDat                    element data
   */
  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) = 0;

  /**
   * @Brief:         
   * 
   * @param outputName 
   * @param outMatrix 
   */
  void AppendNodalOutput(const std::string&outputName, const Matrix&outMatrix);

  /**
   * @Brief:
   * 
   */
  void CommitHistory();

  virtual void GetMassMatrix(std::shared_ptr<ElementData>&elemDat) {}

protected:
  /**
   * @Brief: Compute Number of Element Dof 
   * 
   * @return int 
   */
  inline int DofCount(){
    return m_nodes.size() * m_dofType.size();
  }

protected:
  double m_rho = 0.; 
  std::vector<std::string> m_dofType;                                 // Element Dof Type
  std::shared_ptr<MaterialManager> m_mat;                             // Mateirals
  nlohmann::json m_props;                                             // Whole modele Properties
  std::vector<int> m_nodes;                                           // Element Node Index
  std::unordered_map<std::string, std::vector<double>> m_history;     // History Data
  std::unordered_map<std::string, std::vector<double>> m_current;     // Current Data

protected:
  std::vector<std::vector<double>> xi;
  std::vector<double> weight;
  int order = 0;                                         // the Order of Gauss Integration
  std::string method = "Gauss";                          // the Method of Integration
  Matrix outputData;
};

#endif // ELEMENT_H
