#ifndef TRUSS_H
#define TRUSS_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class Truss : public Element
{
public:
  /**
   * @Brief: Construct a new Truss object
   * 
   * @param elemNode 
   * @param modelProps 
   */
  Truss(const std::vector<int> &elemNode, const nlohmann::json &modelProps);
  
  /**
   * @Brief: Destroy the Truss object
   * 
   */
  ~Truss();

  /**
   * @Brief: Get the Tangent Stiffness object
   * 
   * @param elemDat 
   */
  virtual void GetTangentStiffness(std::shared_ptr<ElementData>&elemDat) override;

private:
  void GetStrain(double &epsilon, double &dEpsilon,
                 const std::vector<double> &a, const std::vector<double> &a0);

  std::vector<double> GetBMatrix(const std::vector<double> &a);

  std::vector<std::vector<double>> GetNonLinearStiffMatrix(const double &sigma,
                                                           const double &area);

private:
  double m_l0 = 0.;
  double m_E = 0.;
  double m_area = 0.;
};

ReflectRegister(Truss, const std::vector<int> &, const nlohmann::json &)

#endif //TRUSS_H