#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <nlohmann/json.hpp>
#include <materials/MaterialManager.h>
#include <util/Math.h>

using namespace Eigen;

/**
 * @Brief:  element data data structure
 * 
 */

class ElementShapeFunctions;

struct ElementData
{
  inline ElementData(VectorXd&elemState, VectorXd&elemDstate){
    m_state = elemState;
    m_Dstate = elemDstate;
    
    m_fint.setZero(elemState.size());
    m_lumped.setZero(elemState.size());
    m_acce.setZero(elemState.size());
    m_velo.setZero(elemState.size());

    m_stiff.setZero(elemState.size(), elemState.size());
    m_mass.setZero(elemState.size(), elemState.size());
  }

  std::vector<std::string> m_outLabel;           // output variable name
  VectorXd m_state;                       // element state variable vector such as displacement
  VectorXd m_Dstate;                      // increment of element state variable vector
  VectorXd m_velo;                        // 
  VectorXd m_acce;
  VectorXd m_fint;                        // element internal force vector
  VectorXd m_lumped;                      // 
  MatrixXd m_mass;                        // element mass matrix
  MatrixXd m_stiff;                       // element stiffness matrix
  MatrixXd m_coords;                      // element node coordinates
  MatrixXd m_outputData;                  // output variable data
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

  inline VectorXd GetHistoryParameter(const std::string&name){
    return m_history[name];
  }
  
  inline void SetHistoryParameter(const std::string&name, const VectorXd &value){
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
  void AppendNodalOutput(const std::string&outputName, const MatrixXd&outMatrix);

  /**
   * @Brief:
   * 
   */
  void CommitHistory();

  /**
   * @Brief: Compute Mass Matrix
   * 
   * @param elemDat 
   */
  virtual void GetMassMatrix(std::shared_ptr<ElementData>&elemDat);

  /**
   * @Brief: Compute the time step size
   * 
   * @param dtK1 
   * @param elemDistortion 
   */
  inline void SetTimeIncrementPara(double dtK1 = 0., double elemDistortion = 0.){
    m_dtK1 = dtK1;
    m_elemDistortion = elemDistortion;
  }

  inline void ReturnTimeIncrementPara(double &dtK1, double &elemDistortion){
    dtK1 = m_dtK1;
    elemDistortion = m_elemDistortion;
  }

protected:
  /**
   * @Brief: Compute Number of Element Dof 
   * 
   * @return int 
   */
  inline int DofCount(){
    return m_nodes.size() * m_dofType.size();
  }

  void GetNMatrix(const Eigen::VectorXd &h);

  virtual void ComputeElemTimeStep(const std::shared_ptr<ElementShapeFunctions> &res,
                                   const MatrixXd &elemNodeCoords,
                                   const VectorXd &elemNodeDisp,
                                   const double detJac);

  virtual void HourGlassTech(std::shared_ptr<ElementData>&elemDat,
                             const VectorXd &elemNodeDisp,
                             const std::shared_ptr<ElementShapeFunctions> &res,
                             const MatrixXd &pHpX);

protected:
  bool m_reductedIntegration = false;
  double m_rho = 0.;
  double m_waveSpeed = 0.;
  double m_dtK1 = 1.e6;
  double m_elemDistortion = 1.;
  double m_vol = 0;
  std::vector<std::string> m_dofType;                                 // Element Dof Type
  std::shared_ptr<MaterialManager> m_mat;                             // Mateirals
  nlohmann::json m_props;                                             // Whole modele Properties
  std::vector<int> m_nodes;                                           // Element Node Index
  std::unordered_map<std::string, VectorXd> m_history;                // History Data
  std::unordered_map<std::string, VectorXd> m_current;                // Current Data

protected:
  MatrixXd xi;
  VectorXd weight;
  int order = 0;                                                      // the Order of Gauss Integration
  std::string method = "Gauss";                                       // the Method of Integration
  MatrixXd outputData;

private:
  MatrixXd N;
};

#endif // ELEMENT_H
