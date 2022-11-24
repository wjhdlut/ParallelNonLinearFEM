#ifndef ELEMENTSET_H
#define ELEMENTSET_H

#include <memory>

#include <fem/NodeSet.h>
#include <nlohmann/json.hpp>
#include <elements/Element.h>

class ElementSet
{
public:
  ElementSet(std::shared_ptr<NodeSet> &nodes, const nlohmann::json &props);
  ~ElementSet();

public:
  void ReadFromFile(const std::string&fileName);

private:
  void Add(const int elemId, const std::string &modelName, const std::vector<int> &elementNodes);

private:
  std::shared_ptr<NodeSet> m_nodes;
  nlohmann::json m_props;
  std::unordered_map<int, std::shared_ptr<Element>> m_elem;
  std::unordered_map<std::string, std::vector<int>> m_groups;
};

#endif // ELEMENTSET_H