#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>
#include <vector>

#include <nlohmann/json.hpp>
#include <materials/MaterialManager.h>

class Element
{
public:
  Element(const std::vector<int> &elemNodes, const nlohmann::json &modelProps);
  virtual ~Element();

  inline std::vector<int> GetNodes(){
    return m_nodes;
  }

private:

protected:
  std::vector<std::string> m_dofType;
  std::shared_ptr<MaterialManager> m_mat;
  nlohmann::json m_props;
  std::vector<int> m_nodes;
};

#endif // ELEMENT_H
