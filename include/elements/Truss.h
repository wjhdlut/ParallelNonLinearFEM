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
   * @param elemShape
   * @param elemNode 
   * @param modelProps 
   */
  Truss(const std::string &elemShape,
        const std::vector<int> &elemNode,
        const nlohmann::json &modelProps);
  
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

  /**
   * @Brief: Compute the Mass Matrix
   * 
   * @param elemDat 
   */
  virtual void GetMassMatrix(std::shared_ptr<ElementData> &elemDat) override{}

private:
  /**
   * @Brief: Compute the Strain Vector
   * 
   * @param epsilon 
   * @param dEpsilon 
   * @param a 
   * @param a0 
   */
  void GetStrain(double &epsilon, double &dEpsilon,
                 const VectorXd &a, const VectorXd &a0);

  /**
   * @Brief: Compute the Strain Matrix B
   * 
   * @param a 
   */
  void GetBMatrix(const VectorXd &a);

  /**
   * @Brief: Initialize Some Variables
   * 
   */
  void Initialize();

  /**
   * @Brief: Compute the NonLinear Stiffnes Matrix
   * 
   * @param sigma 
   * @param area 
   * @return MatrixXd 
   */
  MatrixXd GetNonLinearStiffMatrix(const double &sigma, const double &area);

private:
  double m_l0 = 0.;
  double m_E = 0.;
  double m_area = 0.;

  double epsilon = 0.;
  double dEpsilon = 0.;
  double dSigma = 0.;

  VectorXd B;
  VectorXd a;
  VectorXd Da;
  VectorXd a0;
  VectorXd sigma;
  MatrixXd KL;
  MatrixXd KNL;
};

ReflectRegister(Truss, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif //TRUSS_H