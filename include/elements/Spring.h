/**
 * @File Name:     Spring.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

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
  Spring(const std::string &elemShape,
         const std::vector<int> &elemNode,
         const nlohmann::json &modelProps);

  /**
   * @Brief: Destroy the Spring object
   * 
   */
  ~Spring();

  /**
   * @Brief: Compute the Tangent Stiffness
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
  double m_k = 0.;

private:
  double elong = 0.;
  double Fs = 0.;

  VectorXd a;
  VectorXd Da;
};

ReflectRegister(Spring, const std::string &,
                const std::vector<int> &, const nlohmann::json &)

#endif // SPRING_H