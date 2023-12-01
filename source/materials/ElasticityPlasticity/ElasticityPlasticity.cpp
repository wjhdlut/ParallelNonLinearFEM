/**
 * @File Name:     ElasticityPlasticity.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <iostream>
#include <util/Tools.h>
#include <materials/ElasticityPlasticity/ElasticityPlasticity.h>
#include <util/ObjectFactory.h>

ElasticityPlasticity::ElasticityPlasticity(const nlohmann::json &matProps)
                     : BaseMaterial(matProps)
{
  Initialize();
}

ElasticityPlasticity::~ElasticityPlasticity()
{
}

void ElasticityPlasticity::ReadYieldStressStrainData()
{
  // std::cout << std::setw(4) << m_props << std::endl;
  if(!m_props.contains("yieldStressStrain")) return;

  std::vector<std::string> dataPair = m_props.at("yieldStressStrain");
  for(int i = 0; i < dataPair.size(); i = i + 2)
  {
    m_yieldStress.emplace_back(std::make_pair<double, double>
                            (std::stod(dataPair[i]),std::stod(dataPair[i+1])));
  }
}

void ElasticityPlasticity::ComputeDMatrix()
{
  if(m_yieldFlag)
  {

  }
  else
    m_D = m_lineMat->GetTangMatrix();
}

VectorXd ElasticityPlasticity::GetStress(const std::shared_ptr<Kinematics> &kin,
                                         const VectorXd &increDisp,
                                         const MatrixXd &dphi)
{
  VectorXd stress = m_lineMat->GetStress(kin);
  if("Jaumann" == m_rateType)
    StressRotation(stress, increDisp, dphi);
  
  if(m_yieldRule == nullptr){
    throw " Please Assign the Yile Type";
  }
  double J2 = m_yieldRule->CompYieldFunction(stress);
  double yieldStress = GetYieldStress(m_plasticStrain);

  // Plastic Yield
  if((J2 - yieldStress)/yieldStress > m_tol)
  {
    // Plastic Step: Apply return mapping - use Newton-Raphson algorithm
    //               to solve the return mapping equation
    double denom = 0., dPlasticStrain = 0.;
    while(true)
    {
      denom = -3. * m_G - GetHardModuli(m_plasticStrain);
      dPlasticStrain += (J2 - yieldStress)/denom;

      m_plasticStrain += (J2 - yieldStress)/denom;
      yieldStress = GetYieldStress(m_plasticStrain);

      if((J2 - 3. * m_G * dPlasticStrain - yieldStress)/yieldStress < m_tol)
        break;
    }
    // Update Stress
    //  hydrostatic pressure
    double hydPre = m_oneVec.dot(stress) / 3.;
    // stress deviator
    VectorXd devStress = stress - m_oneVec * hydPre;
    double tempFactor = 1. - 3. * m_G * dPlasticStrain / J2;
    stress = tempFactor * devStress + m_oneVec * hydPre;
  }

  return stress;
}

void ElasticityPlasticity::StressRotation(VectorXd &stress,
                                          const VectorXd &increDisp,
                                          const MatrixXd &dphi)
{
  MatrixXd dispMatrix = Math::ConvertVecToMat(dphi.rows(), dphi.cols(), increDisp);
  
  /**----------------------------------------------
   *  the matrix form of velocity gradient
   *    veloGradient = [pu/px, pu/py, pu/pz 
   *                    pv/px, pv/py, pv/pz
   *                    pw/px, pw/py, pw/pz]
   * ---------------------------------------------- */
  MatrixXd veloGradient =  dphi * dispMatrix;;
  
  // compute dissymmetric velocity gradient
  VectorXd rotStress;
  if(veloGradient.size() == 3){
    Vector3d rotVec;
    rotVec(0) = (veloGradient(0, 1) - veloGradient(1, 0))/2.;
    rotVec(1) = (veloGradient(1, 2) - veloGradient(2, 1))/2.;
    rotVec(2) = (veloGradient(0, 2) - veloGradient(2, 0))/2.;
  
    rotStress = VectorXd::Zero(6);
    rotStress(0) = 2. * (+rotVec(0) * stress(3) + rotVec(2) * stress(5));
    rotStress(1) = 2. * (-rotVec(0) * stress(3) + rotVec(1) * stress(4));
    rotStress(2) = 2. * (-rotVec(2) * stress(5) - rotVec(3) * stress(4));

    rotStress(3) = rotVec(0) * (stress(1) - stress(0)) + rotVec(1) * stress(5) + rotVec(3) * stress(4);
    rotStress(4) = rotVec(2) * (stress(2) - stress(1)) - rotVec(0) * stress(5) - rotVec(2) * stress(3);
    rotStress(5) = rotVec(3) * (stress(2) - stress(0)) + rotVec(0) * stress(4) - rotVec(2) * stress(3);
  }
  stress += rotStress;
}

void ElasticityPlasticity::Initialize()
{
  m_lineMat = std::make_shared<LinearElasticity>(m_props);
  m_E   = SetMaterialParamter("E");
  m_nu  = SetMaterialParamter("nu");
  m_rho = SetMaterialParamter("rho");
  m_plaMod = SetMaterialParamter("tanMod");
  
  m_waveSpeed = sqrt(m_E*(1.-m_nu)/(1+m_nu)*(1-2.*m_nu)*m_rho);

  if(m_props.contains("rateType"))
    m_rateType = m_props.at("rateType");

  if(m_props.contains("yieldFunction")){
    std::string yieldFunction = m_props.at("yieldFunctions");
    m_yieldRule = ObjectFactory::CreateObject<YieldRule>(yieldFunction);
  }

  m_plasticStrain = 0.;

  // Analyse Type
  if(m_props.contains("analyseType"))
  {
    if("PlaneStrain" == m_props.at("analyseType"))
    {
      m_oneVec = VectorXd::Zero(3);
      m_oneVec(0) = 1., m_oneVec(1) = 1.;

      m_oneMat = MatrixXd::Zero(3, 3);
      m_oneMat(0, 0) = 1., m_oneMat(1, 1) = 1., m_oneMat(2, 2) = 0.5;

      m_devMat = MatrixXd::Zero(3, 3);
      m_devMat = m_oneMat - m_oneVec * m_oneVec.transpose() / 3.;
    }
    if("PlaneStress" == m_props.at("analyseType"))
    {
      // To Do
    }
  }
  else
  {
    m_oneVec = VectorXd::Zero(6);
    m_oneVec(0) = 1., m_oneVec(1) = 1., m_oneVec(2) = 1.;

    m_oneMat = MatrixXd::Zero(6, 6);
    m_oneMat(0, 0) = 1.0, m_oneMat(1, 1) = 1.0, m_oneMat(2, 2) = 1.0;
    m_oneMat(3, 3) = 0.5, m_oneMat(4, 4) = 0.5, m_oneMat(5, 5) = 0.5;
    
    m_devMat = MatrixXd::Zero(6, 6);
    m_devMat = m_oneMat - m_oneVec * m_oneVec.transpose() / 3.;
  }

  ReadYieldStressStrainData();
}

double ElasticityPlasticity::GetYieldStress(const double plasticStrain)
{
  int index = 0;
  for(index = 0; index < m_yieldStress.size(); index++)
  {
    if(plasticStrain > m_yieldStress.at(index).first){
      continue;
    }
    else{
      /* x <= x0  --------> f(x) = f(x0);
         this case is impossible for Plastic Problem
      */
      if(0 == index)
        return m_yieldStress.at(index).second;
      else{
        /*  y(i) - y(i-1)      y - y(i-1)
         * --------------- = -------------
         *  x(i) 0- x(i-1)     x - x(i-1)
         */
        double k1 = (m_yieldStress.at(index).second-m_yieldStress.at(index-1).second);
        double k2 = (m_yieldStress.at(index).first-m_yieldStress.at(index-1).first);
        double k3 = (plasticStrain - m_yieldStress.at(index - 1).first);
        return m_yieldStress.at(index-1).second + k1 / k2 * k3;
      }
    }
  }
  return m_yieldStress.at(index).second;
}

double ElasticityPlasticity::GetHardModuli(const double plasticStrain)
{
  int index = 0;
  for(index = 0; index < m_yieldStress.size(); index++)
  {
    if(plasticStrain > m_yieldStress.at(index).first)
      continue;
    else{
      if(0 == index)
        return 0.;
      else{
        double k1 = (m_yieldStress.at(index).second-m_yieldStress.at(index-1).second);
        double k2 = (m_yieldStress.at(index).first-m_yieldStress.at(index-1).first);
        return k1/k2;
      }
    }
  }
  return 0.;
}