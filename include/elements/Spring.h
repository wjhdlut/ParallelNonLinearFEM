#ifndef SPRING_H
#define SPRING_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class Spring : public Element
{
public:
  /**
   * @Brief: Construct a new Spring object
   * 
   * @param elemNode 
   * @param modelProps 
   */
  Spring(const std::vector<int> &elemNode, const nlohmann::json &modelProps);

  /**
   * @Brief: Destroy the Spring object
   * 
   */
  ~Spring();

  /**
   * @Brief: Get the Tangent Stiffness object
   * 
   * @param elemDat 
   */
  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

  virtual void GetMassMatrix(std::shared_ptr<ElementData> &elemDat) override{}

private:
  double m_k = 0.;

private:
  double elong = 0.;
  double Fs = 0.;

  VectorXd a;
  VectorXd Da;
};

ReflectRegister(Spring, const std::vector<int> &, const nlohmann::json &)

#endif // SPRING_H