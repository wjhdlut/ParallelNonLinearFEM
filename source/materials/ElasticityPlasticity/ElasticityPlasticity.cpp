#include <materials/ElasticityPlasticity/ElasticityPlasticity.h>

ElasticityPlasticity::ElasticityPlasticity(const nlohmann::json &matProps) : BaseMaterial(matProps)
{
  Initialize();
}

ElasticityPlasticity::~ElasticityPlasticity()
{
  
}

void ElasticityPlasticity::ComputeDMatrix()
{
  m_D = m_lineMat->GetTangMatrix();
}

VectorXd ElasticityPlasticity::GetStress(const std::shared_ptr<Kinematics> &kin,
                                         const VectorXd &increDisp,
                                         const MatrixXd &dphi)
{
  VectorXd stress = m_lineMat->GetStress(kin);
  if("Jaumann" == m_rateType)
    StressRotation(stress, increDisp, dphi);

  double meanStress = m_oneVec.dot(stress);
  VectorXd devStress = stress - meanStress * m_oneVec;

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
    ratio * devStress;
  }

  return devStress + meanStress * m_oneVec;
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
  m_yieldStress = SetMaterialParamter("yieldStress");
  m_plaMod = SetMaterialParamter("tanMod");
  
  m_waveSpeed = sqrt(m_E*(1.-m_nu)/(1+m_nu)*(1-2.*m_nu)*m_rho);

  if(m_props.contains("rateType"))
    m_rateType = m_props.at("rateType");

  m_oneVec << 1., 1., 1., 0., 0., 0.;
}