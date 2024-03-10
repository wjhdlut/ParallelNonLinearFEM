/**
 * @File Name:     EquivantStrainBase.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2024-01-23
 * 
 * @Copyright Copyright (c) 2024 JianHuaWang
 * 
 */

#include "../include/EquivantStrainBase.h"

EquivantStrainBase::EquivantStrainBase(const nlohmann::json&props) : m_props(props)
{
  Initialize();
}

EquivantStrainBase::~EquivantStrainBase()
{}

void EquivantStrainBase::Initialize()
{
  if(m_props.contains("analyseType"))
  {
    if("PlaneStrain"== m_props.at("analyseType"))
      m_planeStrainFlag = true;
    if("PlaneStress" == m_props.at("analyseType"))
      m_planeStressFlag = true;
  }
}