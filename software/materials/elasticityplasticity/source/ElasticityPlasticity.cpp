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

#include "../include/ElasticityPlasticity.h"
#include "../../../util/include/Tools.h"
#include "../../../util/include/ObjectFactory.h"

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
  m_D = m_lineMat->GetTangMatrix();
  
  if(m_yieldFlag)
  {
    // double qTrial = m_yieldRule->CompYieldFunction(stress);
    double qTrial = m_J2;
    // double qTrial = J2;
    // double H = GetHardModuli(m_accumPlasticStrain);
    double devStress3 = 0.; 
    if(m_planeStrainFlag)
      devStress3 = -(m_devStress[0] + m_devStress[1]);

    double normDevStress = 0.;
    if(3 == m_oneVec.size())
      normDevStress = sqrt(pow(m_devStress(0), 2) + pow(m_devStress(1), 2)
                           + pow(devStress3, 2.) + 2. * (pow(m_devStress(2), 2)));
    if(6 == m_oneVec.size()){
      normDevStress = sqrt(pow(m_devStress(0), 2) + pow(m_devStress(1), 2) + pow(m_devStress(2), 2)
                    + 2. * (pow(m_devStress(3), 2) + pow(m_devStress(4), 2) + pow(m_devStress(5), 2)));
    }

    // Compute the Consistent Tangent Modulus
    // std::cout << "m_devMat = " << m_devMat << std::endl;
    // std::cout << "m_devStress = \n" << m_devStress << std::endl;
    // std::cout << "A = " << m_devStress * m_devStress.transpose() << std::endl;
    m_D = m_D - 2. * m_G * (3. * m_G * m_dPlasticMultiplier) / qTrial * m_devMat 
        + 6. * m_G * m_G * (m_dPlasticMultiplier / qTrial - 1. / (3. * m_G + m_H))/ (normDevStress * normDevStress)
        * m_devStress * m_devStress.transpose();
    
  }
  // std::cout << "m_D = \n" << m_D << std::endl;
}

VectorXd ElasticityPlasticity::GetStress(const std::shared_ptr<Kinematics> &kin)
{
  m_yieldFlag     = false;
  // Compute Trial Stress Based Displacement
  VectorXd stressTrial = CompElasticTrialStress(kin);
  // std::cout << "stressTrial = " << stressTrial.transpose() << std::endl;
  // std::cout << "strain = " << kin->strain.transpose() << std::endl;
  // std::cout << "incrementalStrain = " << kin->incremStrain.transpose() << std::endl;
  // if("Jaumann" == m_rateType)
    // StressRotation(stress, increDisp, dphi);
  
  GetHistoryParameter(m_accumPlasticStrain, "AccumPlasticStrain");
  double yieldStress = GetYieldStress(m_accumPlasticStrain);
  double J2 = m_yieldRule->CompYieldFunction(stressTrial);
  m_J2 = J2;

  // Plastic Yield
  if((J2 - yieldStress)/yieldStress > m_tol)
  {
    m_yieldFlag     = true;
    m_updateDMatrix = true;
    // Plastic Step: Apply return mapping - use Newton-Raphson algorithm
    //               to solve the return mapping equation
    double denom = 0.;
    m_dPlasticMultiplier = 0.;
    while(true)
    {
      m_H = GetHardModuli(m_accumPlasticStrain);
      denom = 3. * m_G + m_H;
      m_dPlasticMultiplier += (J2 - yieldStress)/denom;

      m_accumPlasticStrain += (J2 - yieldStress)/denom;
      yieldStress = GetYieldStress(m_accumPlasticStrain);

      if((J2 - 3. * m_G * m_dPlasticMultiplier - yieldStress)/yieldStress < m_tol)
        break;
    }
    // Update Stress
    // Hydrostatic Stress
    double stress3 = 0.;
    if(m_planeStrainFlag)
      stress3 = m_nu * (stressTrial[0] + stressTrial[1]);
    m_hydPre = (m_oneVec.dot(stressTrial) + stress3) / 3.;
    
    // Deviatoric Stress Component
    VectorXd devStressTrial = stressTrial - m_oneVec * m_hydPre;
    double tempFactor = 1. - 3. * m_G * m_dPlasticMultiplier / J2;

    m_devStress = tempFactor * devStressTrial;
    
    stressTrial = m_devStress + m_hydPre * m_oneVec;    
  }
  // std::cout << "stress = " << stressTrial.transpose() << std::endl;
  
  SetHistoryParameter("AccumPlasticStrain", m_accumPlasticStrain);
  SetHistoryParameter("stress", stressTrial);

  return stressTrial;
}

VectorXd ElasticityPlasticity::CompElasticTrialStress(const std::shared_ptr<Kinematics> &kin)
{
  VectorXd stress;

  GetHistoryParameter(stress, "stress");

  // std::cout << "incremStrain = " << kin->incremStrain.transpose() << std::endl;

  return stress + m_lineMat->GetTangMatrix() * kin->incremStrain;
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
  // Read Basic Material Parameters
  m_lineMat = std::make_shared<LinearElasticity>(m_props);
  m_E       = SetMaterialParamter("E");
  m_nu      = SetMaterialParamter("nu");
  m_rho     = SetMaterialParamter("rho");
  m_plaMod  = SetMaterialParamter("tanMod");
  m_G       = m_E /(2. * (1+m_nu));
  
  m_waveSpeed = sqrt(m_E*(1.-m_nu)/(1+m_nu)*(1-2.*m_nu)*m_rho);

  if(m_props.contains("rateType"))
    m_rateType = m_props.at("rateType");

  if(m_props.contains("yieldFunction")){
    std::string yieldFunction = m_props.at("yieldFunction");
    m_yieldRule = ObjectFactory::CreateObject<YieldRule>(yieldFunction, m_props);
  }
  
  if(m_yieldRule == nullptr){
    std::cout << "Catch Exception: "
              << " Please Assign the Yile Type for the Elasto-Plastic Material!!!"
              << std::endl;
    exit(-1);
  }

  // m_accumPlasticStrain = 0.;
  // Initialize the Effective Plastic Strain
  SetHistoryParameter("AccumPlasticStrain", 0.);

  // Analyse Type
  if(m_planeStrainFlag)
  {
    m_devStress = VectorXd::Zero(3);
    m_oneVec    = VectorXd::Zero(3);
    m_oneVec(0) = 1., m_oneVec(1) = 1.;

    m_oneMat = MatrixXd::Zero(3, 3);
    m_oneMat(0, 0) = 1., m_oneMat(1, 1) = 1., m_oneMat(2, 2) = 0.5;

    m_devMat = MatrixXd::Zero(3, 3);
    m_devMat = m_oneMat - m_oneVec * m_oneVec.transpose() / 3.;

    SetHistoryParameter("stress", VectorXd::Zero(3));
  }
  else if(m_planeStressFlag)
  {
    // To Do
    m_devStress = VectorXd::Zero(3);
    SetHistoryParameter("stress", VectorXd::Zero(3));
  }
  else
  {
    // For 3D Case
    m_devStress = VectorXd::Zero(6);
    m_oneVec = VectorXd::Zero(6);
    m_oneVec(0) = 1., m_oneVec(1) = 1., m_oneVec(2) = 1.;

    m_oneMat = MatrixXd::Zero(6, 6);
    m_oneMat(0, 0) = 1.0, m_oneMat(1, 1) = 1.0, m_oneMat(2, 2) = 1.0;
    m_oneMat(3, 3) = 0.5, m_oneMat(4, 4) = 0.5, m_oneMat(5, 5) = 0.5;
    
    m_devMat = MatrixXd::Zero(6, 6);
    m_devMat = m_oneMat - m_oneVec * m_oneVec.transpose() / 3.;

    SetHistoryParameter("stress", VectorXd::Zero(6));
  }

  ReadYieldStressStrainData();
}

double ElasticityPlasticity::GetYieldStress(const double accumPlasticStrain)
{
  int index = 0;
  for(index = 0; index < m_yieldStress.size(); index++)
  {
    if(accumPlasticStrain > m_yieldStress.at(index).first){
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
        double k3 = (accumPlasticStrain - m_yieldStress.at(index - 1).first);
        return m_yieldStress.at(index-1).second + k1 / k2 * k3;
      }
    }
  }
  return m_yieldStress.at(index).second;
}

double ElasticityPlasticity::GetHardModuli(const double accumPlasticStrain)
{
  int index = 0;
  for(index = 0; index < m_yieldStress.size(); index++)
  {
    if(accumPlasticStrain > m_yieldStress.at(index).first)
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