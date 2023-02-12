
#ifndef DOFSPACE_H
#define DOFSPACE_H

#include <map>
#include <fem/ElementSet.h>
#include <petscksp.h>

class DofSpace
{
public:
  DofSpace(std::shared_ptr<ElementSet> elements, std::shared_ptr<NodeSet> nodes);
  ~DofSpace();

  void ReadFromFile(const std::string&fileName);

  inline void SetConstrainFactor(double fac){
    m_constrainedFac = fac;
  }

  inline std::vector<std::string> GetDofType(){
    return m_dofTypes;
  }

  inline int GetForType(const int nodeIds, const std::string&dofType){
    int indexRow = std::distance(m_IDmap.begin(), std::find(m_IDmap.begin(), m_IDmap.end(), nodeIds));
    int indexLine = std::distance(m_dofTypes.begin(), std::find(m_dofTypes.begin(), m_dofTypes.end(), dofType));
    return m_dofs[indexRow][indexLine];
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
   * @Brief: solve equation when K is scale
   * 
   * @param K 
   * @param df 
   * @param da 
   */
  PetscErrorCode Solve(double K, Vec&df, Vec&da);

  /**
   * @Brief: Compute the Norm of dF Vector for the value at non-Constrained Location 
   * 
   * @param r                   [in]  dF Vector
   * @param error               [out] the Norm Value
   * @return PetscErrorCode 
   */
  PetscErrorCode Norm(Vec &r, double &error);

public:
  std::vector<std::vector<int>> m_dofs;

private:
  void Constrain(const int&nodeId, const std::string&dofType, const double&value);

  PetscErrorCode GetConstraintsMatrix(Mat&C);

private:
  std::vector<std::string> m_dofTypes;
  std::map<int, std::vector<double>> *m_nodeCoords;
  std::unordered_map<int, double> m_constrained;
  double m_constrainedFac = 1.;
  std::vector<int> m_IDmap;
};

#endif // DOFSPACE_H
