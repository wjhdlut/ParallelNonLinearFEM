/**
 * @File Name:     ElementSet.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2024-01-11
 * 
 * @Copyright Copyright (c) 2024 JianHuaWang
 * 
 */

#ifndef ELEMENTSET_H
#define ELEMENTSET_H

#include <memory>
#include <petscmat.h>
#include "NodeSet.h"
#include "../../nlohmann/json.hpp"
#include "../../elements/include/Element.h"

class Element;

class ElementSet
{
public:
  /**
   * @Brief: Construct a new Element Set object
   * 
   * @param nodes 
   * @param props 
   */
  ElementSet(std::shared_ptr<NodeSet> &nodes, const nlohmann::json &props);
  
  /**
   * @Brief: Destroy the Element Set object
   * 
   */
  ~ElementSet();

public:
  /**
   * @Brief: Read Element Data From File
   * 
   * @param fileName 
   */
  void ReadFromFile(const std::string&fileName);

  /**
   * @Brief: Get the Dof Type of Element
   * 
   * @return std::vector<std::string> 
   */
  std::vector<std::string> GetDofType();

  /**
   * @Brief: Compute and Assemble Stiffness Matrix
   * 
   * @param A 
   * @param B 
   */
  void AssembleTangentStiffness(Mat &A, Vec &B);
  
  /**
   * @Brief: Compute and Assemble Internal Force
   * 
   * @param B 
   */
  void AssembleInternalForce(Vec &B);
  
  /**
   * @Brief: Compute and Assemble Internal Force
   * 
   * @param A 
   * @param B 
   * @return PetscErrorCode 
   */
  PetscErrorCode AssembleMassMatrix(Mat &A, Vec &B);

  void CommitHistory();

  inline void GetCurrentTimeIncrementPara(double dtK1,
                                          double elemDistortion = 0)
  {
    m_dtK1 = dtK1;
    m_elemDistortion = elemDistortion;
  }

  inline void ReturnCurentTimeIncrementPara(double &dtK1,
                                            double &elemDistortion)
  {
    dtK1 = m_dtK1;
    elemDistortion = m_elemDistortion;
  }

  /**
   * @Brief:  return number of elements in group with name groupName
   * 
   * @param groupName 
   * @return int 
   */
  int ElemGroupCount(const std::string&groupName);

  /**
   * @Brief: return element id vector in group with name groupName
   * 
   * @param elemGroupName 
   * @return std::vector<int> 
   */
  std::vector<int> IterElementGroup(const std::string&elemGroupName);

  inline std::unordered_map<int, std::shared_ptr<Element>> GetElementPtr(){
    return m_elem;
  }

private:
  void Add(const int elemId, const std::string &modelName,
           const std::string &elemShape, const std::vector<int> &elementNodes);
  
  /**
   * @Brief:  Assemble the element Matrix (such as Element Stiffness Matrix
   *  and Mass Matrix) into the Global Matrix
   * 
   * @param A 
   * @param B 
   * @param rank 
   * @param action 
   * @return PetscErrorCode 
   */
  PetscErrorCode AssembleMatrix(Mat &A, Vec&B, const int rank, const std::string&action);

  inline void ResetDtime(){
    m_dtK1 = 1.e6;
  }

private:
  std::shared_ptr<NodeSet> m_nodes;
  nlohmann::json m_props;
  std::unordered_map<int, std::shared_ptr<Element>> m_elem;
  std::unordered_map<std::string, std::vector<int>> m_groups; // [groupName, elementIndex]

private:
  std::shared_ptr<Element> elemPtr;
  std::shared_ptr<ElementData> elemData;
  std::vector<int> elemNodes;
  std::vector<int> elemDofs;
  VectorXd elemState;
  VectorXd elemDstate;
  MatrixXd elemCoords;

private:
  double m_dtK1 = 1.e6;
  double m_elemDistortion = 1.;
};

#endif // ELEMENTSET_H