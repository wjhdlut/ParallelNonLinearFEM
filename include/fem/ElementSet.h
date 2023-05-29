#ifndef ELEMENTSET_H
#define ELEMENTSET_H

#include <memory>

#include <fem/NodeSet.h>
#include <nlohmann/json.hpp>
#include <elements/Element.h>
#include <petscmat.h>

class Element;

class ElementSet
{
public:
  ElementSet(std::shared_ptr<NodeSet> &nodes, const nlohmann::json &props);
  ~ElementSet();

public:
  void ReadFromFile(const std::string&fileName);

  std::vector<std::string> GetDofType();

  void AssembleTangentStiffness(Mat &A, Vec &B);
  
  void AssembleInternalForce(Vec &B);
  
  void AssembleMassMatrix(Mat &A, Vec &B);

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
  void Add(const int elemId, const std::string &modelName, const std::vector<int> &elementNodes);
  
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
  std::vector<double> elemState;
  std::vector<double> elemDstate;
  std::vector<std::vector<double>> elemCoords;

private:
  double m_dtK1 = 1.e6;
  double m_elemDistortion = 1.;
};

#endif // ELEMENTSET_H