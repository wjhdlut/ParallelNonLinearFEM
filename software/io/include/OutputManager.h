/**
 * @File Name:     OutputManager.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-12-04
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <vector>

#include "../../util/include/BaseModule.h"

class OutputManager
{
public:
  /**
   * @Brief: Construct a new Output Manager object
   * 
   */
  OutputManager();

  /**
   * @Brief: Destroy the Output Manager object
   * 
   */
  ~OutputManager();

  /**
   * @Brief: output related variables
   * 
   */
  void Run();

private:
  std::vector<std::shared_ptr<BaseModule>> m_outman;
};

#endif // OUTPUTMANAGER_H