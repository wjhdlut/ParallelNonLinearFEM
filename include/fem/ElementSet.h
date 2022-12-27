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
  void AssembleTangentStiffness(Mat &A, Vec&B);
  void AssembleInternalForce(Vec &B);
  void AssembleMassMatrix(Mat &A, Vec&B);

private:
  void Add(const int elemId, const std::string &modelName, const std::vector<int> &elementNodes);
  
  PetscErrorCode AssembleMatrix(Mat &A, Vec&B, const int rank, const std::string&action);

private:
  std::shared_ptr<NodeSet> m_nodes;
  nlohmann::json m_props;
  std::unordered_map<int, std::shared_ptr<Element>> m_elem;
  std::unordered_map<std::string, std::vector<int>> m_groups;
};

#endif // ELEMENTSET_H