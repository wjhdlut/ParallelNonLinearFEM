/**
 * @File Name:     Element.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "../../nlohmann/json.hpp"
#include "../../materials/include/MaterialManager.h"
#include "../../util/include/Math.h"

using namespace Eigen;

/**
 * @Brief:  element data data structure
 * 
 */

class ElementShapeFunctions;


/**
 * @Brief: Element Basic Data
 * 
 */
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

  VectorXd m_state;                                                       // element state variable vector such as displacement
  VectorXd m_Dstate;                                                      // increment of element state variable vector
  VectorXd m_velo;                                                        // Velocity
  VectorXd m_acce;                                                        // Acceleration
  VectorXd m_fint;                                                        // Element internal force vector
  VectorXd m_lumped;                                                      // Mass Matrix
  MatrixXd m_mass;                                                        // Element mass matrix
  MatrixXd m_stiff;                                                       // Element stiffness matrix
  MatrixXd m_coords;                                                      // Element node coordinates
  std::unordered_map<std::string, MatrixXd> m_outputData;                 // Output variable data
};

/**
 * @Brief: father class of element
 * 
 */
class Element
{
public:
  /**
   * @Brief: Construct a new Element Class
   * 
   * @param elemNodes 
   * @param modelProps 
   */
  Element(const std::vector<int> &elemNodes,
          const nlohmann::json &modelProps);
  
  /**
   * @Brief: Destroy the Element Class
   * 
   */
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

  /**
   * @Brief: Return Element Time Increment Parameter
   * 
   * @param dtK1 
   * @param elemDistortion 
   */
  inline void ReturnTimeIncrementPara(double &dtK1, double &elemDistortion){
    dtK1 = m_dtK1;
    elemDistortion = m_elemDistortion;
  }
  
  /**
   * @Brief: Compute Components of the Equivalent Nodal Loads based on Edge Load Data
   * 
   * @param nodeForcePres 
   * @return VectorXd 
   */
  virtual VectorXd CompEquivalentNodeForce(const std::unordered_map<int, VectorXd> &nodeCoord,
                                           const std::unordered_map<int, std::vector<double>> &nodeForcePres);

protected:
  /**
   * @Brief: Compute Number of Element Dof 
   * 
   * @return int 
   */
  inline int DofCount(){
    return m_nodes.size() * m_dofType.size();
  }
  
  /**
   * @Brief: Compute Shape Function Matrix
   * 
   * @param h 
   */
  void GetNMatrix(const Eigen::VectorXd &h);

  /**
   * @Brief: Compute Element Time Step
   * 
   * @param res 
   * @param elemNodeCoords 
   * @param elemNodeDisp 
   * @param detJac 
   */
  virtual void ComputeElemTimeStep(const MatrixXd &elemNodeCoords,
                                   const VectorXd &elemNodeDisp,
                                   const double detJac);
  
  /**
   * @Brief: Compute Hourglass Force For Explicit Dynamic Problem
   * 
   * @param elemDat 
   * @param elemNodeDisp 
   * @param res 
   * @param pHpX 
   */
  virtual void HourGlassTech(std::shared_ptr<ElementData>&elemDat,
                             const VectorXd &elemNodeDisp,
                             const MatrixXd &pHpX);
  
  /**
   * @Brief: Set the Dof Type m_dofType for Different Element
   * 
   */
  void SetDofType(const std::string &elemShape);
 
 /**
  * @Brief: Compute Gauss Point Coordinate
  * 
  * @param elemShape 
  */
  void CompGaussPointCoord(const std::string &elemShape);
  
  /**
   * @Brief: Initialize History Variables
   * 
   */
  void InitializeHistoryVariables();
 
  /**
   * @Brief: Check Whether a Given Set of Local Element Node Correspond to
   *         One of the Element Boundaries (Edges in 2-D and Faces in 3-D).
   *         if it does, Returns the Local Node Numbers Ordered for Numerical
   *         Integration on Boundary.
   * 
   * @param nodeChk 
   * @param nodeForcePres 
   * @return std::vector<int>  [out]  the Local Node Numbers Ordered
   */
  std::vector<int> CheckNodeBoundary(std::vector<int> &nodeChk,
                                     const std::unordered_map<int, std::vector<double>> &nodeForcePres);

  /**
   * @Brief: Insert Element Output Data to Specify Variables
   * 
   * @param elemOutData 
   * @param name 
   * @param outputData 
   */
  void InsertElemOutputData(std::unordered_map<std::string, MatrixXd> &elemOutData,
                            const std::string &name,
                            const MatrixXd &outputData = MatrixXd::Zero(0, 0));
  
  /**
   * @Brief: Compute the Strain Matrix for AxiSymmetric Problem
   * 
   * @param B                    [in/out] Strain Matrix
   * @param h                    [in]    Shape funtion
   * @param elemDat              [in]    Element Data
   */
  double GetBMatrixForAxiSymmetry(MatrixXd &B,
                                  const MatrixXd &h,
                                  const MatrixXd &nodeCoord);

private:
  /**
   * @Brief: Initialize Some Basic Variables
   * 
   * @param modelProps 
   */
  void Initialize(const nlohmann::json &modelProps);

  /**
   * @Brief: Initialize Stress Vector at Gauss Point in Element
   * 
   */
  void InitializeStress();

  /**
   * @Brief: Initialize State Variables Such as Plastic Multiplier
   * 
   */
  void InitializeStateVariable();

  /**
   * @Brief: Set the Analyse Type
   * 
   * @param iter 
   */
  void SetAnalyseType(const nlohmann::json &modelProps);
  
  /**
   * @Brief: Set the Material Parameters
   * 
   * @param modelProps 
   */
  void SetMaterialParameters(const nlohmann::json &modelProps);

protected:
  bool   m_reductedIntegration = false;
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
  int order = 0;                                                      // the Order of Gauss Integration
  std::string m_analyseType;                                          // the Analyse Type
  std::string m_axiSymmetry;                                          // the Axis of Symmetry for the Axisymetric Problem
  MatrixXd xi;
  VectorXd weight;
  std::string method = "Gauss";                                       // the Method of Integration
  std::shared_ptr<ElementShapeFunctions> m_elemShapePtr;              // the Element Shape Pointer
  MatrixXd outputData;                                                // the Element Output Data 
  std::unordered_map<int, std::vector<int>> m_nodeOrdered;            // the Element Node Ordered

private:
  MatrixXd N;
  const double Pi = 3.1415926;
};

#endif // ELEMENT_H