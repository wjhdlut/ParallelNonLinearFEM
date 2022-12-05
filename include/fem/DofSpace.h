
#ifndef DOFSPACE_H
#define DOFSPACE_H

#include <map>
#include <fem/ElementSet.h>

class DofSpace
{
public:
  DofSpace(std::shared_ptr<ElementSet> elements, std::shared_ptr<NodeSet> nodes);
  ~DofSpace();

  void ReadFromFile(const std::string&fileName);

  inline void SetConstrainFactor(double fac){
    m_constrainedFac = fac;
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

public:
  std::vector<std::vector<int>> m_dofs;

private:
  void Constrain(const int&nodeId, const std::string&dofType, const double&value);

private:
  std::vector<std::string> m_dofTypes;
  std::map<int, std::vector<double>> *m_nodeCoords;
  std::unordered_map<int, double> m_constrained;
  double m_constrainedFac = 1.;
  std::vector<int> m_IDmap;
};

#endif // DOFSPACE_H
