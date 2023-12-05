/**
 * @File Name:     OutputWriter.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef OUTPUTWRITER_H
#define OUTPUTWRITER_H

#include "../../util/include/BaseModule.h"
#include "../../util/include/ObjectFactory.h"

class OutputWriter : public BaseModule
{
public:
  /**
   * @Brief: Construct a new Output Writer object
   * 
   * @param props 
   */
  OutputWriter(const nlohmann::json&props);

  /**
   * @Brief: Destroy the Output Writer object
   * 
   */
  ~OutputWriter();

  /**
   * @Brief:
   * 
   */
  virtual void Run() override;

private:
  bool m_onScreen = false;
};

ReflectRegister(OutputWriter, const nlohmann::json&)

#endif // OUTPUTWRITER_H