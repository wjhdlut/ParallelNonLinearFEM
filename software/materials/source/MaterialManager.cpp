/**
 * @File Name:     MaterialManager.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <iostream>

#include "../include/MaterialManager.h"
#include "../../util/include/ObjectFactory.h"

MaterialManager::MaterialManager(const nlohmann::json &matProps)
{
  std::string matType = matProps.at("type");
  m_mat = ObjectFactory::CreateObject<BaseMaterial>(matType, matProps);
  if(nullptr == m_mat){
    std::cout <<  matType << " Material Model Created Failed!!" << std::endl;
    exit(-1);
  }
  iIter = -1;
}

MaterialManager::~MaterialManager()
{}

void MaterialManager::Reset()
{
  iIter = -1;
}

VectorXd MaterialManager::GetStress(const std::shared_ptr<Kinematics>&kin,
                                    int iSam)
{
  if(-1 == iSam){
    iIter += 1;
    iSam = iIter;
  }

  m_mat->SetIter(iSam);

  return m_mat->GetStress(kin);
}