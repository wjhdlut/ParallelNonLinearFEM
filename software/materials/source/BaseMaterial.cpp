/**
 * @File Name:     BaseMaterial.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include "../include/BaseMaterial.h"

BaseMaterial::BaseMaterial(const nlohmann::json &props) : m_props(props)
{
}

BaseMaterial::~BaseMaterial()
{}

void BaseMaterial::SetIter(int iIter)
{
  m_iIter = iIter;
}

int BaseMaterial::SetHistoryParameter(const std::string &name, double value)
{
  if(-1 == m_iIter) 
  {
    m_initHistory[name] = value;
    return 0;
  }

  if(m_current.size() == m_iIter) m_current.emplace_back(m_initHistory);

  m_current[m_iIter].at(name) = value;

  return 0;
}

double BaseMaterial::GetHistoryParameter(const std::string&name)
{
  if(0 == m_history.size()) return m_initHistory[name];
  else return m_history[m_iIter].at(name);
}

void BaseMaterial::CommitHistory()
{
  m_history.clear();
  for(auto h : m_current)
    m_history.emplace_back(h);
}

VectorXd BaseMaterial::GetStress(const std::shared_ptr<Kinematics> &kin,
                                 const VectorXd &increDisp, const MatrixXd &dhpi){
  return m_D * kin->strain;
}

double BaseMaterial::SetMaterialParamter(const std::string &name)
{
  if(!m_props.contains(name)) return 0.;
  if(m_props.at(name).is_string()){
    std::string E = m_props.at(name);
    return std::stod(E);
  }
  else{
    return m_props.at(name);
  }
}

void BaseMaterial::TangentDMatrixToCurrent(const MatrixXd &F)
{
  MatrixXd T;
  GetTransMatrix(T, F);
  m_D = T * m_D * T.transpose();
}


void BaseMaterial::TangentDMatrixToInitial(const MatrixXd &F)
{
  MatrixXd T;
  GetTransMatrix(T, F);
  m_D = T.inverse() * m_D * T.inverse().transpose();
}

void BaseMaterial::GetTransMatrix(MatrixXd &T, const MatrixXd &F)
{
  if(2 == F.rows())
  {
    T.resize(3, 3);
  }
  if(3 == F.rows())
  {
    T.resize(6, 6);

    T(0, 0) = F(0, 0) * F(0, 0);
    T(0, 1) = F(0, 1) * F(0, 1);
    T(0, 2) = F(0, 2) * F(0, 2);
    T(0, 3) = F(0, 0) * F(0, 1) + F(0, 1) * F(0, 0);
    T(0, 4) = F(0, 1) * F(0, 2) + F(0, 2) * F(0, 1);
    T(0, 5) = F(0, 0) * F(0, 2) + F(0, 2) * F(0, 1);

    T(1, 0) = F(1, 0) * F(1, 0);
    T(1, 1) = F(1, 1) * F(1, 1);
    T(1, 2) = F(1, 2) * F(1, 2);
    T(1, 3) = F(1, 0) * F(1, 1) + F(1, 1) * F(1, 0);
    T(1, 4) = F(1, 1) * F(1, 2) + F(1, 2) * F(1, 1);
    T(1, 5) = F(1, 0) * F(1, 2) + F(1, 2) * F(1, 0);

    T(2, 0) = F(2, 0) * F(2, 0);
    T(2, 1) = F(2, 1) * F(2, 1);
    T(2, 2) = F(2, 2) * F(2, 2);
    T(2, 3) = F(2, 0) * F(2, 1) + F(2, 1) * F(2, 0);
    T(2, 4) = F(2, 1) * F(2, 2) + F(2, 2) * F(2, 1);
    T(2, 5) = F(2, 0) * F(2, 2) + F(2, 2) * F(2, 0);

    T(3, 0) = F(0, 0) * F(1, 0);
    T(3, 1) = F(0, 1) * F(1, 1);
    T(3, 2) = F(0, 2) * F(1, 2);
    T(3, 3) = F(0, 0) * F(1, 1) + F(0, 1) * F(1, 0);
    T(3, 4) = F(0, 1) * F(1, 2) + F(0, 2) * F(1, 1);
    T(3, 5) = F(0, 0) * F(1, 2) + F(0, 2) * F(1, 0);

    T(4, 0) = F(1, 0) * F(2, 0);
    T(4, 1) = F(1, 1) * F(2, 1);
    T(4, 2) = F(1, 2) * F(2, 2);
    T(4, 3) = F(1, 0) * F(2, 1) + F(1, 1) * F(2, 0);
    T(4, 4) = F(1, 1) * F(2, 2) + F(1, 2) * F(2, 1);
    T(4, 5) = F(1, 0) * F(2, 2) + F(1, 2) * F(2, 0);

    T(5, 0) = F(0, 0) * F(2, 0);
    T(5, 1) = F(0, 1) * F(2, 1);
    T(5, 2) = F(0, 2) * F(2, 2);
    T(5, 3) = F(0, 0) * F(2, 1) + F(0, 1) * F(2, 0);
    T(5, 4) = F(0, 1) * F(2, 2) + F(0, 2) * F(2, 1);
    T(5, 5) = F(0, 0) * F(2, 2) + F(0, 2) * F(2, 0);
  }
}