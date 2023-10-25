/**
 * @File Name:     Hexa8ShapeFunctions.cpp
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

#include <elements/shapefunctions/Hexa8ShapeFunctions.h>

#include <iostream>

Hexa8ShapeFunctions::Hexa8ShapeFunctions()
{
  Initialize();
};

Hexa8ShapeFunctions::~Hexa8ShapeFunctions()
{}

void Hexa8ShapeFunctions::Initialize()
{
  H     = VectorXd::Zero(8);
  pHpxi = MatrixXd::Zero(8, 3);

  numOfStress = 6;
  
  m_face = MatrixXi::Zero(6, 4);
  m_face(0, 0) = 1, m_face(0, 1) = 2, m_face(0, 2) = 3, m_face(0, 3) = 4;
  m_face(1, 0) = 5, m_face(1, 1) = 6, m_face(1, 2) = 7, m_face(1, 3) = 8;
  m_face(2, 0) = 1, m_face(2, 1) = 2, m_face(2, 2) = 6, m_face(2, 3) = 5;
  m_face(3, 0) = 2, m_face(3, 1) = 3, m_face(3, 2) = 7, m_face(3, 3) = 6;
  m_face(4, 0) = 3, m_face(4, 1) = 4, m_face(4, 2) = 8, m_face(4, 3) = 7;
  m_face(5, 0) = 4, m_face(5, 1) = 1, m_face(5, 2) = 5, m_face(5, 3) = 8;
  
  m_h(0, 0) = 1., m_h(0, 1) =-1, m_h(0, 2) = 1, m_h(0, 3) = -1;
  m_h(1, 0) = 1., m_h(1, 1) = 1, m_h(1, 2) =-1, m_h(1, 3) = -1;
  m_h(2, 0) = 1., m_h(2, 1) =-1, m_h(2, 2) =-1, m_h(2, 3) =  1;
  m_h(3, 0) = 1., m_h(3, 1) =-1, m_h(3, 2) = 1, m_h(3, 3) = -1;

  m_ss(0, 0) = 2., m_ss(0, 1) =-2., m_ss(0, 2) =2., m_ss(0, 3) =-2.;
  m_ss(1, 0) =-2., m_ss(1, 1) =-2., m_ss(1, 2) =2., m_ss(1, 3) = 2.;
  m_ss(2, 0) =-2., m_ss(2, 1) = 2., m_ss(2, 2) =2., m_ss(2, 3) =-2.;
  m_ss(3, 0) = 0., m_ss(3, 1) = 0., m_ss(3, 2) =0., m_ss(3, 3) = 0.;
  
  m_dofType = {"u", "v", "w"};
}

void Hexa8ShapeFunctions::GetShapeFunction(const VectorXd &xi)
{
  if(3 != xi.size()) throw "The isoparamatric coordinate should be 3D for Hexa8 element.";
  
  // compute shape funtion values at gauss point
  H(0) = 0.125*(1.0-xi(0))*(1.0-xi(1))*(1.0-xi(2));
  H(1) = 0.125*(1.0+xi(0))*(1.0-xi(1))*(1.0-xi(2));
  H(2) = 0.125*(1.0+xi(0))*(1.0+xi(1))*(1.0-xi(2));
  H(3) = 0.125*(1.0-xi(0))*(1.0+xi(1))*(1.0-xi(2));
  H(4) = 0.125*(1.0-xi(0))*(1.0-xi(1))*(1.0+xi(2));
  H(5) = 0.125*(1.0+xi(0))*(1.0-xi(1))*(1.0+xi(2));
  H(6) = 0.125*(1.0+xi(0))*(1.0+xi(1))*(1.0+xi(2));
  H(7) = 0.125*(1.0-xi(0))*(1.0+xi(1))*(1.0+xi(2));

  // Calculate derivatives of shape functions 
  pHpxi(0, 0) = -0.125*(1.0-xi(1))*(1.0-xi(2));
  pHpxi(1, 0) =  0.125*(1.0-xi(1))*(1.0-xi(2));
  pHpxi(2, 0) =  0.125*(1.0+xi(1))*(1.0-xi(2));
  pHpxi(3, 0) = -0.125*(1.0+xi(1))*(1.0-xi(2));
  pHpxi(4, 0) = -0.125*(1.0-xi(1))*(1.0+xi(2));
  pHpxi(5, 0) =  0.125*(1.0-xi(1))*(1.0+xi(2));
  pHpxi(6, 0) =  0.125*(1.0+xi(1))*(1.0+xi(2));
  pHpxi(7, 0) = -0.125*(1.0+xi(1))*(1.0+xi(2));

  pHpxi(0, 1) = -0.125*(1.0-xi(0))*(1.0-xi(2));
  pHpxi(1, 1) = -0.125*(1.0+xi(0))*(1.0-xi(2));
  pHpxi(2, 1) =  0.125*(1.0+xi(0))*(1.0-xi(2));
  pHpxi(3, 1) =  0.125*(1.0-xi(0))*(1.0-xi(2));
  pHpxi(4, 1) = -0.125*(1.0-xi(0))*(1.0+xi(2));
  pHpxi(5, 1) = -0.125*(1.0+xi(0))*(1.0+xi(2));
  pHpxi(6, 1) =  0.125*(1.0+xi(0))*(1.0+xi(2));
  pHpxi(7, 1) =  0.125*(1.0-xi(0))*(1.0+xi(2));

  pHpxi(0, 2) = -0.125*(1.0-xi(0))*(1.0-xi(1));
  pHpxi(1, 2) = -0.125*(1.0+xi(0))*(1.0-xi(1));
  pHpxi(2, 2) = -0.125*(1.0+xi(0))*(1.0+xi(1));
  pHpxi(3, 2) = -0.125*(1.0-xi(0))*(1.0+xi(1));
  pHpxi(4, 2) =  0.125*(1.0-xi(0))*(1.0-xi(1));
  pHpxi(5, 2) =  0.125*(1.0+xi(0))*(1.0-xi(1));
  pHpxi(6, 2) =  0.125*(1.0+xi(0))*(1.0+xi(1));
  pHpxi(7, 2) =  0.125*(1.0-xi(0))*(1.0+xi(1));
}

VectorXd Hexa8ShapeFunctions::HourGlassTech(std::shared_ptr<ElementData> &elemDat,
                                            const VectorXd &elemNodeDisp,
                                            const double &c,
                                            const nlohmann::json &hourGlassPara,
                                            const MatrixXd &pHpX)
{
  double para = hourGlassPara.at("para");
  if("STANDARD" == hourGlassPara.at("type"))
    return HourGlassStand(elemDat, c, para);
  else if("Flanagan-Belytschko" == hourGlassPara.at("type"))
    return HourGlassFlangan(elemDat, elemNodeDisp, c, pHpX);

  return VectorXd::Zero(0);
}

VectorXd Hexa8ShapeFunctions::HourGlassFlangan(std::shared_ptr<ElementData> &elemDat,
                                               const VectorXd &elemNodeDisp,
                                               const double &c,
                                               const MatrixXd &pHpX)
{
  double x3478, x2358, x1467, x1256;
  beat.resize(3, 4);
  MatrixXd nodeDispMat = Map<const MatrixXd>(elemNodeDisp.data(), 3, 8).transpose();
  for(int i = 0; i < 3; i++){
    x3478 = (elemDat->m_coords(2, i) + nodeDispMat(2, i))
          - (elemDat->m_coords(3, i) + nodeDispMat(3, i))
          - (elemDat->m_coords(6, i) + nodeDispMat(6, i))
          + (elemDat->m_coords(7, i) + nodeDispMat(7, i));
    
    x2358 = (elemDat->m_coords(1, i) + nodeDispMat(1, i))
          - (elemDat->m_coords(2, i) + nodeDispMat(2, i))
          - (elemDat->m_coords(4, i) + nodeDispMat(4, i))
          + (elemDat->m_coords(7, i) + nodeDispMat(7, i));

    x1467 = (elemDat->m_coords(0, i) + nodeDispMat(0, i))
          - (elemDat->m_coords(3, i) + nodeDispMat(3, i))
          - (elemDat->m_coords(5, i) + nodeDispMat(5, i))
          + (elemDat->m_coords(6, i) + nodeDispMat(6, i));

    x1256 = (elemDat->m_coords(0, i) + nodeDispMat(0, i))
          - (elemDat->m_coords(1, i) + nodeDispMat(1, i))
          - (elemDat->m_coords(4, i) + nodeDispMat(4, i))
          + (elemDat->m_coords(5, i) + nodeDispMat(5, i));

    // x3478 = elemDat->m_coords(2, i) - elemDat->m_coords(3, i)
    //       - elemDat->m_coords(6, i) + elemDat->m_coords(7, i);
    // x2358 = elemDat->m_coords(1, i) - elemDat->m_coords(2, i)
    //       - elemDat->m_coords(4, i) + elemDat->m_coords(7, i);
    // x1467 = elemDat->m_coords(0, i) - elemDat->m_coords(3, i)
    //       - elemDat->m_coords(5, i) + elemDat->m_coords(6, i);
    // x1256 = elemDat->m_coords(0, i) - elemDat->m_coords(1, i)
    //       - elemDat->m_coords(4, i) + elemDat->m_coords(5, i);

    beat(i, 0) = x1467 - x2358;
    beat(i, 1) = x1467 + x2358;
    beat(i, 2) = x1256 - x3478;
    beat(i, 3) = x1256 + x3478;
  }
  
  MatrixXd gama = MatrixXd::Zero(4, 8);
  for(int j = 0; j < 4; j++)
    for(int k = 0; k < 4; k++)
    {
      gama(j, k) = m_h(j, k);
      for(int i = 0; i < 3; i++)
        gama(j, k) = gama(j, k) - beat(i, j) * pHpX(k, i);
    }
  
  for(int j = 0; j < 4; j++){
    gama(j, 4) = m_ss(j, 0) - gama(j, 2);
    gama(j, 5) = m_ss(j, 1) - gama(j, 3);
    gama(j, 6) = m_ss(j, 2) - gama(j, 0);
    gama(j, 7) = m_ss(j, 3) - gama(j, 1);
  }

  int index = 0;
  MatrixXd g = MatrixXd::Zero(pHpX.cols(), gama.rows());
  for(int i = 0; i < pHpX.cols(); i++){
    for(int k = 0; k < gama.rows(); k++){
      for(int j = 0; j < pHpX.rows(); j++){
        index = j * pHpX.cols() + i;
        g(i, k) += elemDat->m_velo[index] * gama(k, j);
      }
    }
  }
  
  VectorXd temp = VectorXd::Zero(elemDat->m_fint.size());
  MatrixXd tempMat = g * gama;
  for(int i = 0; i < pHpX.rows(); i++){
    for(int j = 0; j < pHpX.cols(); j++){
      index = pHpX.cols() * i + j;
      temp(index) -= c * tempMat(j, i);
    }
  }

  return temp;
}

VectorXd Hexa8ShapeFunctions::HourGlassStand(std::shared_ptr<ElementData> &elemDat,
                                             const double &c,
                                             const double &para)
{
  double v3478, v2358, v1467, v1256;
  hgr.resize(3, 4);
  for(int i = 0; i < 3; i++){
    v3478 = elemDat->m_velo( 6+i) - elemDat->m_velo( 9+i)
          - elemDat->m_velo(18+i) + elemDat->m_velo(21+i);
    v2358 = elemDat->m_velo( 3+i) - elemDat->m_velo( 6+i)
          - elemDat->m_velo(12+i) + elemDat->m_velo(21+i);
    v1467 = elemDat->m_velo( 0+i) - elemDat->m_velo( 9+i)
          - elemDat->m_velo(15+i) + elemDat->m_velo(18+i);
    v1256 = elemDat->m_velo( 0+i) - elemDat->m_velo( 3+i)
          - elemDat->m_velo(12+i) + elemDat->m_velo(15+i);

    hgr(i, 0) = v1467 - v2358;
    hgr(i, 1) = v1467 + v2358;
    hgr(i, 2) = v1256 - v3478;
    hgr(i, 3) = v1256 + v3478;
  }

  double ah = c, al = 100. * para;
  for(int i = 0; i < 3; i++){
    for(int k = 0; k < 4; k++){
      hgr(i, k) = hgr(i, k) * (ah + abs(al*hgr(i, k)));
    }
  }
  double hap, ham, hbp, hbm;
  VectorXd temp(elemDat->m_fint.size());
  for(int i = 0; i < 3; i++){
    hap = hgr(i, 0) + hgr(i, 3);
    ham = hgr(i, 0) - hgr(i, 3);
    hbp = hgr(i, 1) + hgr(i, 2);
    hbm = hgr(i, 1) - hgr(i, 2);

    temp( 0+i) = -hap - hbp;
    temp( 3+i) =  hap - hbm;
    temp( 6+i) = -hap + hbp;
		temp( 9+i) =  hap + hbm;
		temp(12+i) = -ham + hbp;
		temp(15+i) =  ham + hbm;
		temp(18+i) = -ham - hbp;
		temp(21+i) =  ham - hbm;
  }

  return temp;
}

double Hexa8ShapeFunctions::ComputeElemTimeStep(double &dtK1, double &elemDistortion,
                                                const MatrixXd &elemNodeCoords,
                                                const VectorXd &elemNodeDisp,
                                                const double detJac,  const double waveSpeed)
{
  int k1, k2, k3, k4;
  double areal = 1.0e20, aream = 0.;
  double e, g, f, atest;
  double x13, x24, fs, ft;
  for(int iFace = 0; iFace < m_face.rows(); iFace++)
  {
      k1 = m_face(iFace, 0) - 1;
      k2 = m_face(iFace, 1) - 1;
      k3 = m_face(iFace, 2) - 1;
      k4 = m_face(iFace, 3) - 1;
    
    e = 0., f = 0., g = 0.;
    for(int iDof = 0; iDof < m_dofType.size(); iDof++){
      x13 = (elemNodeCoords(k3, iDof) + elemNodeDisp(m_dofType.size() * k3 + iDof))
          - (elemNodeCoords(k1, iDof) + elemNodeDisp(m_dofType.size() * k1 + iDof));
      x24 = (elemNodeCoords(k4, iDof) + elemNodeDisp(m_dofType.size() * k4 + iDof))
          - (elemNodeCoords(k2, iDof) + elemNodeDisp(m_dofType.size() * k2 + iDof));
      
      // x13 = nodeCoords(k3, iDof) - nodeCoords(k1, iDof);
      // x24 = nodeCoords(k4, iDof) - nodeCoords(k2, iDof);

      fs = x13 - x24;
      ft = x13 + x24;
      
      e += fs * fs;
      f += fs * ft;
      g += ft * ft;
    }
    atest = e*g - f*f;

    aream = std::max(atest, aream);
    areal = std::min(atest, areal);
  }
  
  // double vol = 8. * detJac;
  double at = areal / aream;
  double dt = 4 * 8. * detJac / sqrt(aream);
  
  dt /= waveSpeed;
  
  dtK1 = std::min(dtK1, dt);
  elemDistortion = std::min(elemDistortion, at);

  return 8. * detJac;
}

std::unordered_map<int, std::vector<int>> Hexa8ShapeFunctions::SetElemNodeOrdered()
{
  std::unordered_map<int, std::vector<int>> elemNodeOrdered;
  elemNodeOrdered.insert(std::pair<int, std::vector<int>>(1, {1, 2, 3, 0, 0, 0, 0, 0}));
  elemNodeOrdered.insert(std::pair<int, std::vector<int>>(2, {0, 0, 1, 2, 3, 0, 0, 0}));
  elemNodeOrdered.insert(std::pair<int, std::vector<int>>(3, {0, 0, 0, 0, 1, 2, 3, 0}));
  elemNodeOrdered.insert(std::pair<int, std::vector<int>>(4, {3, 0, 0, 0, 0, 0, 1, 2}));

  return elemNodeOrdered;
}