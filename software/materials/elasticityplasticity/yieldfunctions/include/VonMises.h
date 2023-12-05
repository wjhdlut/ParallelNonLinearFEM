/**
 * @File Name:     VonMises.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-26
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "YieldRule.h"
#include "../../../../util/include/ObjectFactory.h"

class VonMises : public YieldRule
{
public:
  /**
   * @Brief: Construct a new Von Mises Material object
   * 
   */
  VonMises(const nlohmann::json &matProps);

  /**
   * @Brief: Destroy the Von Mises Material object
   * 
   */
  ~VonMises();
  
  /**
   * @Brief:         
   * 
   * @param stress 
   */
  virtual double CompYieldFunction(const VectorXd &stress) override;

private:
  VectorXd m_oneVec;
};
ReflectRegister(VonMises, const nlohmann::json &)