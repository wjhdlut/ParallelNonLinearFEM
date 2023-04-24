#include <materials/ElasticityPlasticity/ElasticityPlasticity.h>

ElasticityPlasticity::ElasticityPlasticity(const nlohmann::json &matProps) : BaseMaterial(matProps)
{
  lineMat = std::make_shared<LinearElasticity>(matProps);
  m_E   = SetMaterialParamter("E");
  m_nu  = SetMaterialParamter("nu");
  m_rho = SetMaterialParamter("rho");
  m_yieldStress = SetMaterialParamter("yieldStress");
  m_plaMod = SetMaterialParamter("tanMod");
  m_waveSpeed = sqrt(m_E*(1.-m_nu)/(1+m_nu)*(1-2.*m_nu)*m_rho);

  if(matProps.contains("rateType"))
    m_rateType = matProps.at("rateType");
}

ElasticityPlasticity::~ElasticityPlasticity()
{
  
}

void ElasticityPlasticity::ComputeDMatrix()
{
  m_D = lineMat->GetTangMatrix();
}

std::vector<double> ElasticityPlasticity::GetStress(const std::shared_ptr<Kinematics> &kin,
                                                    const std::vector<double> &increDisp,
                                                    const Matrix &dphi)
{
  std::vector<double> stress = lineMat->GetStress(kin);
  if("Jaumann" == m_rateType)
    StressRotation(stress, increDisp, dphi);

  double meanStress = Math::VecDot(m_oneVec, stress);
  std::vector<double> devStress = Math::VecAdd(-1., stress, Math::VecScale(meanStress, m_oneVec));

  double J2 = 0.5 * (devStress[0] * devStress[0] + devStress[1] * devStress[1] + devStress[2] * devStress[2]) + 
                     devStress[3] * devStress[3] + devStress[4] * devStress[4] + devStress[5] * devStress[5];
  
  J2 = sqrt(3. * J2);

  if(J2 > m_yieldStress)
  {
    double depeff = (J2 - m_yieldStress) / (3. * m_G + m_plaMod);
    double epeff_ = epeff_ + depeff;

    m_yieldStress = m_yieldStress + m_plaMod*depeff;
    
    double ratio = m_yieldStress / J2;
    
    J2 = J2 * ratio;
    Math::VecScale(ratio, devStress);
  }

  return Math::VecAdd(1., devStress, Math::VecScale(meanStress, m_oneVec));
}

void ElasticityPlasticity::StressRotation(std::vector<double> &stress,
                                          const std::vector<double> &increDisp,
                                          const Matrix &dphi)
{
  Matrix dispMatrix = Math::VecReshape(dphi.size(), dphi[0].size(), increDisp);
  
  /**----------------------------------------------
   *  the matrix form of velocity gradient
   *    veloGradient = [pu/px, pu/py, pu/pz 
   *                    pv/px, pv/py, pv/pz
   *                    pw/px, pw/py, pw/pz]
   * ---------------------------------------------- */
  Matrix veloGradient =  Math::MatrixATransMultB(dphi, dispMatrix);
  
  // compute dissymmetric velocity gradient
  std::vector<double> rotStress;
  if(veloGradient.size() == 3){
    std::vector<double> rotVec = {(veloGradient[0][1] - veloGradient[1][0])/2.,
                                  (veloGradient[1][2] - veloGradient[2][1])/2.,
                                  (veloGradient[0][2] - veloGradient[2][0])/2.};
  
    rotStress.resize(6, 0.);
    rotStress[0] = 2. * (+rotVec[0] * stress[3] + rotVec[2] * stress[5]);
    rotStress[1] = 2. * (-rotVec[0] * stress[3] + rotVec[1] * stress[4]);
    rotStress[2] = 2. * (-rotVec[2] * stress[5] - rotVec[3] * stress[4]);

    rotStress[3] = rotVec[0] * (stress[1] - stress[0]) + rotVec[1] * stress[5] + rotVec[3] * stress[4];
    rotStress[4] = rotVec[2] * (stress[2] - stress[1]) - rotVec[0] * stress[5] - rotVec[2] * stress[3];
    rotStress[5] = rotVec[3] * (stress[2] - stress[0]) + rotVec[0] * stress[4] - rotVec[2] * stress[3];
  }
  stress = Math::VecAdd(1., stress, rotStress);
}