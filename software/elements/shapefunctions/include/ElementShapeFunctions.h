/**
 * @File Name:     ElementShapeFunctions.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef ELEMENTSHAPEFUNCTIONS_H
#define ELEMENTSHAPEFUNCTIONS_H

#include <vector>
#include <string>
#include "../../../nlohmann/json.hpp"
#include "../../include/Element.h"

/**
 * @Brief:   Base Class to Compute Element Shape Functions
 * 
 */

/* pHpxi = [pH1pxi1 pH1pxi2 pH1pxi3,
            pH2pxi1 pH2pxi2 pH2pxi3，
            ...,
            pHnpxi1 pHnpxi2 pHnpxi3]*/

struct ElementData;

class ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct a new Element Shape Functions object
   * 
   */
  ElementShapeFunctions(){}

  /**
   * @Brief: Destroy the Element Shape Functions object
   * 
   */
  virtual ~ElementShapeFunctions(){}

  /**
   * @Brief: Compute the Shape Function Virtual Method
   * 
   * @param xi                [in]  Gauss Coordinate
   */
  virtual void GetShapeFunction(const VectorXd &xi) = 0;
  
  /**
   * @Brief: Compute the Boundary Shape Function
   * 
   * @param boundaryH 
   * @param pboundaryHpxi 
   */
  inline virtual void GetBoundaryShapeFunction(VectorXd &boundaryH,
                                               MatrixXd &pboundaryHpxi,
                                               const VectorXd &boundaryXi){

  }

  /**
   * @Brief: Get the Boundary Integration Point
   * 
   * @param boundaryXi 
   * @param boundaryWeight 
   */
  inline virtual void GetBoundaryIntegrationPoint(MatrixXd &boundaryXi, VectorXd &boundaryWeight){

  }
  
  /**
   * @Brief: Compute Element Hourglass Force
   * 
   * @param elemDat 
   * @param elemNodeDisp 
   * @param c 
   * @param para 
   * @param pHpX 
   * @return VectorXd 
   */
  virtual VectorXd HourGlassTech(std::shared_ptr<ElementData> &elemDat,
                                 const VectorXd &elemNodeDisp,
                                 const double &c,
                                 const nlohmann::json &para,
                                 const MatrixXd &pHpX)
  {
    return VectorXd::Zero(0);
  }
  
  /**
   * @Brief: Return Element Matrix
   * 
   * @return MatrixXi* 
   */
  inline MatrixXi *ReturnFaceMatrix(){
    return &m_face;
  }
  
  /**
   * @Brief: Return Dof Type m_dofType
   * 
   * @return std::vector<std::string> 
   */
  inline std::vector<std::string> ReturnDofType(){
    return m_dofType;
  }
  
  /**
   * @Brief: Compute Element Time Step
   * 
   * @param dtK1 
   * @param elemDistortion 
   * @param elemNodeCoords 
   * @param elemNodeDisp 
   * @param detJac 
   * @param waveSpeed 
   * @return double 
   */
  virtual double ComputeElemTimeStep(double &dtK1, double &elemDistortion,
                                     const MatrixXd &elemNodeCoords,
                                     const VectorXd &elemNodeDisp,
                                     const double detJac, const double waveSpeed)
  {
    return 0.;
  }

  inline std::unordered_map<int, std::vector<int>> ReturnElemNodeOrdered(){
    return SetElemNodeOrdered();
  }

protected:
  /**
   * @Brief: Inilize Element Variables
   * 
   */
  virtual void Initialize() = 0;
  
  /**
   * @Brief: Set the Element Node Ordered object
   * 
   * @return std::unordered_map<int, std::vector<int>> 
   */
  inline virtual std::unordered_map<int, std::vector<int>> SetElemNodeOrdered(){
    return std::unordered_map<int, std::vector<int>> {};
  }

public:
  int numOfStress;                              // The number of stress
  VectorXd H;                                   // The shape function
  MatrixXd pHpxi;                               // The derivative of shape function about local coordinate system

protected:
  MatrixXi m_face;
  std::vector<std::string> m_dofType;           // Element Dof Type
};

#endif // ELEMENTSHAPEFUNCTIONS_H
