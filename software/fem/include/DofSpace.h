/**
 * @File Name:     DofSpace.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2024-03-14
 * 
 * @Copyright Copyright (c) 2024 JianHuaWang
 * 
 */

#ifndef DOFSPACE_H
#define DOFSPACE_H

#include <map>
#include <petscksp.h>
#include "ElementSet.h"

struct RigidWall{
  int direction = 0;
  double coord = 0.;
};

class DofSpace
{
public:
  /**
   * @Brief: Construct a new DofSpace object
   * 
   * @param elements 
   * @param nodes 
   */
  DofSpace(std::shared_ptr<ElementSet> elements, std::shared_ptr<NodeSet> nodes);

  /**
   * @Brief: Destroy the DofSpace object
   * 
   */
  ~DofSpace();

  /**
   * @Brief: Read Basic Information From a File
   * 
   * @param fileName 
   */
  void ReadFromFile(const std::string&fileName);

  
  inline void SetConstrainFactor(double fac){
    m_constrainedFac = fac;
  }
  
  /**
   * @Brief: Get the Dof Type of Element
   * 
   * @return std::vector<std::string> 
   */
  inline std::vector<std::string> GetDofType(){
    return m_dofTypes;
  }

  inline int GetForType(const int nodeIds, const std::string&dofType){
    int indexRow = std::distance(m_IDmap.begin(), std::find(m_IDmap.begin(), m_IDmap.end(), nodeIds));
    int indexLine = std::distance(m_dofTypes.begin(), std::find(m_dofTypes.begin(), m_dofTypes.end(), dofType));
    return m_dofs[indexRow][indexLine];
  }

  inline std::vector<int> GetForType(const int nodeIds){
    int indexRow = std::distance(m_IDmap.begin(), std::find(m_IDmap.begin(), m_IDmap.end(), nodeIds));
    return m_dofs[indexRow];
  }

  inline std::vector<int> Get(const std::vector<int>&elemNodes){
    std::vector<int> elemDofs;
    for(auto node : elemNodes)
    {
      int indexRow = std::distance(m_IDmap.begin(), std::find(m_IDmap.begin(), m_IDmap.end(), node));
      for(auto dof : m_dofs.at(indexRow))
        elemDofs.emplace_back(dof);
    }
    return elemDofs;
  }

  inline std::vector<int> GetIndex(const std::vector<int>&nodeIDs){
    std::vector<int> index;
    for(auto nodeID : nodeIDs)
      index.emplace_back(GetIndex(nodeID));
    return index;
  }

  inline int GetIndex(int nodeID){
    return std::distance(m_IDmap.begin(), std::find(m_IDmap.begin(), m_IDmap.end(), nodeID));
  }


  /**
   * @Brief:  solve equation sets when K is matrix
   * 
   * @param K  stiffness matrix
   * @param df increment of force vector
   * @param da displacement increment
   */
  PetscErrorCode Solve(Mat&K, Vec&df, Vec&da, KSP&ksp);

  /**
   * @Brief: solve equation when K is Vector
   * 
   * @param K 
   * @param df 
   * @param da 
   */
  PetscErrorCode Solve(const Vec&K, const Vec&df, Vec&da);

  /**
   * @Brief: Compute the Norm of dF Vector for the value at non-Constrained Location 
   * 
   * @param r                   [in]  dF Vector
   * @param error               [out] the Norm Value
   * @return PetscErrorCode 
   */
  PetscErrorCode Norm(const Vec &r, double &error);

  /**
   * @Brief: Compute Some Basic Variables
   * 
   */
  void CompBasicVariable();

  double Converg(const Vec &fext, const Vec &fint);

private:
  void Constrain(const int&nodeId, const std::string&dofType, const double&value);

  /**
   * @Brief: Compute the Constraints Matrix
   * 
   * @param C 
   * @return PetscErrorCode 
   */
  PetscErrorCode GetConstraintsMatrix(Mat&C);

  PetscErrorCode GetConstraintsMatrix();

  /**
   * @Brief: Read Node COnstraint From File
   * 
   * @param fileName 
   */
  void ReadNodeConstraint(const std::string &fileName);

  /**
   * @Brief: Read Rigid Wall Information From FIle
   * 
   * @param fileName 
   */
  void ReadRigidWall(const std::string &fileName);

  /**
   * @Brief: Apply the Rigid Wall COnstraint
   * 
   * @param da 
   */
  void RigidWallConstraint(Vec&da);

  /**
   * @Brief: Transform Matrx (such as the Stiffness or Mass Matrix)
   *         to Cylinderical Coordinate System
   * 
   * @param A 
   */
  PetscErrorCode TransMatToCylinder(Mat&A);

  /**
   * @Brief: Transform the Vector (Internal Force Vector)
   *         to Cylinderical Coordinate System
   * 
   * @param B 
   */
  PetscErrorCode TransVecToCylinder(Vec &B);
  
  /**
   * @Brief: Transform the Vector (Displacement Vector)
   *         to Global Cartesian Coordinate System
   * 
   * @param B 
   * @return PetscErrorCode 
   */
  PetscErrorCode TransVecToGlobalCSY(Vec &B);

  /**
   * @Brief: Initialize Some Basic Information
   * 
   */
  void Initialize();

  /**
   * @Brief: Compute the Transform Matrix
   * 
   */
  PetscErrorCode ComputeTransMatrix();

public:
  std::vector<std::vector<int>> m_dofs;

private:
  std::vector<std::string> m_dofTypes;
  std::map<int, VectorXd>  *m_nodeCoords;
  std::unordered_map<int, double> m_constrained;
  double m_constrainedFac = 1.;
  std::vector<int> m_IDmap;
  std::shared_ptr<RigidWall> m_rigidWall = nullptr;
  Mat m_transMatrix = NULL;
  Mat m_C = NULL;
};

#endif // DOFSPACE_H
