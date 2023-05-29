#include <elements/shapefunctions/Hexa8ShapeFunctions.h>

Hexa8ShapeFunctions::Hexa8ShapeFunctions()
{
  H.resize(8);
  
  std::vector<double> temp(3, 0.);
  pHpxi.resize(8, temp);

  numOfStress = 6;

  m_face = {{1, 2, 3, 4},
            {5, 6, 7, 8},
            {1, 2, 6, 5},
            {2, 3, 7, 6},
            {3, 4, 8, 7},
            {4, 1, 5, 8}};
};

Hexa8ShapeFunctions::~Hexa8ShapeFunctions()
{}

void Hexa8ShapeFunctions::GetShapeFunction(const std::vector<double> &xi)
{
  if(3 != xi.size()) throw "The isoparamatric coordinate should be 3D for Hexa8 element.";
  
  // compute shape funtion values at gauss point
  H[0] = 0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0-xi[2]);
  H[1] = 0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[2]);
  H[2] = 0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[2]);
  H[3] = 0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[2]);
  H[4] = 0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[2]);
  H[5] = 0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0+xi[2]);
  H[6] = 0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0+xi[2]);
  H[7] = 0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[2]);

  // Calculate derivatives of shape functions 
  pHpxi[0][0] = -0.125*(1.0-xi[1])*(1.0-xi[2]);
  pHpxi[1][0] =  0.125*(1.0-xi[1])*(1.0-xi[2]);
  pHpxi[2][0] =  0.125*(1.0+xi[1])*(1.0-xi[2]);
  pHpxi[3][0] = -0.125*(1.0+xi[1])*(1.0-xi[2]);
  pHpxi[4][0] = -0.125*(1.0-xi[1])*(1.0+xi[2]);
  pHpxi[5][0] =  0.125*(1.0-xi[1])*(1.0+xi[2]);
  pHpxi[6][0] =  0.125*(1.0+xi[1])*(1.0+xi[2]);
  pHpxi[7][0] = -0.125*(1.0+xi[1])*(1.0+xi[2]);

  pHpxi[0][1] = -0.125*(1.0-xi[0])*(1.0-xi[2]);
  pHpxi[1][1] = -0.125*(1.0+xi[0])*(1.0-xi[2]);
  pHpxi[2][1] =  0.125*(1.0+xi[0])*(1.0-xi[2]);
  pHpxi[3][1] =  0.125*(1.0-xi[0])*(1.0-xi[2]);
  pHpxi[4][1] = -0.125*(1.0-xi[0])*(1.0+xi[2]);
  pHpxi[5][1] = -0.125*(1.0+xi[0])*(1.0+xi[2]);
  pHpxi[6][1] =  0.125*(1.0+xi[0])*(1.0+xi[2]);
  pHpxi[7][1] =  0.125*(1.0-xi[0])*(1.0+xi[2]);

  pHpxi[0][2] = -0.125*(1.0-xi[0])*(1.0-xi[1]);
  pHpxi[1][2] = -0.125*(1.0+xi[0])*(1.0-xi[1]);
  pHpxi[2][2] = -0.125*(1.0+xi[0])*(1.0+xi[1]);
  pHpxi[3][2] = -0.125*(1.0-xi[0])*(1.0+xi[1]);
  pHpxi[4][2] =  0.125*(1.0-xi[0])*(1.0-xi[1]);
  pHpxi[5][2] =  0.125*(1.0+xi[0])*(1.0-xi[1]);
  pHpxi[6][2] =  0.125*(1.0+xi[0])*(1.0+xi[1]);
  pHpxi[7][2] =  0.125*(1.0-xi[0])*(1.0+xi[1]);
}

std::vector<double> Hexa8ShapeFunctions::HourGlassTech(std::shared_ptr<ElementData> &elemDat,
                                                       const double &c,
                                                       const nlohmann::json &hourGlassPara,
                                                       const std::vector<std::vector<double>> &pHpX)
{
  double para = hourGlassPara.at("para");
  if("STANDARD" == hourGlassPara.at("type"))
    return HourGlassStand(elemDat, c, para);
  else if("Flanagan-Belytschko" == hourGlassPara.at("type"))
    return HourGlassFlangan(elemDat, c, pHpX);

  return std::vector<double>(0.);
}

std::vector<double> Hexa8ShapeFunctions::HourGlassFlangan(std::shared_ptr<ElementData> &elemDat,
                                                          const double &c,
                                                          const std::vector<std::vector<double>> &pHpX)
{
  double x3478, x2358, x1467, x1256;
  std::vector<double> temp(4, 0.);
  beat.resize(3, temp);
  for(int i = 0; i < 3; i++){
    x3478 = elemDat->m_coords[2][i] - elemDat->m_coords[3][i] - elemDat->m_coords[6][i] + elemDat->m_coords[7][i];
    x2358 = elemDat->m_coords[1][i] - elemDat->m_coords[2][i] - elemDat->m_coords[4][i] + elemDat->m_coords[7][i];
    x1467 = elemDat->m_coords[0][i] - elemDat->m_coords[3][i] - elemDat->m_coords[5][i] + elemDat->m_coords[6][i];
    x1256 = elemDat->m_coords[0][i] - elemDat->m_coords[1][i] - elemDat->m_coords[4][i] + elemDat->m_coords[5][i];

    beat[i][0] = x1467 - x2358;
    beat[i][1] = x1467 + x2358;
    beat[i][2] = x1256 - x3478;
    beat[i][3] = x1256 + x3478;
  }
  
  temp.resize(8, 0.);
  std::vector<std::vector<double>> gama(4, temp);
  for(int j = 0; j < 4; j++)
    for(int k = 0; k < 4; k++)
    {
      gama[j][k] = m_h[j][k];
      for(int i = 0; i < 3; i++)
        gama[j][k] = gama[j][k] - beat[i][j] * pHpX[k][i];
    }
  
  for(int j = 0; j < 4; j++){
    gama[j][4] = m_ss[j][0] - gama[j][2];
    gama[j][5] = m_ss[j][1] - gama[j][3];
    gama[j][6] = m_ss[j][2] - gama[j][0];
    gama[j][7] = m_ss[j][3] - gama[j][1];
  }

  int index = 0;
  temp.resize(gama.size(), 0.);
  Matrix g(pHpX[0].size(), temp);
  for(int i = 0; i < pHpX[0].size(); i++){
    for(int k = 0; k < gama.size(); k++){
      for(int j = 0; j < pHpX.size(); j++){
        index = j * pHpX[0].size() + i;
        g[i][k] += elemDat->m_velo[index] * gama[k][j];
      }
    }
  }
  
  temp.resize(elemDat->m_fint.size(), 0.);
  Matrix tempMat = Math::MatrixAMultB(g, gama);
  for(int i = 0; i < pHpX.size(); i++){
    for(int j = 0; j < pHpX[0].size(); j++){
      index = pHpX[0].size() * i + j;
      temp[index] -= c * tempMat[j][i];
    }
  }

  return temp;
}

std::vector<double> Hexa8ShapeFunctions::HourGlassStand(std::shared_ptr<ElementData> &elemDat,
                                                        const double &c,
                                                        const double &para)
{
  double v3478, v2358, v1467, v1256;
  std::vector<double> temp(4, 0.);
  hgr.resize(3, temp);
  for(int i = 0; i < 3; i++){
    v3478 = elemDat->m_velo[6+i] - elemDat->m_velo[9+i] - elemDat->m_velo[18+i] + elemDat->m_velo[21+i];
    v2358 = elemDat->m_velo[3+i] - elemDat->m_velo[6+i] - elemDat->m_velo[12+i] + elemDat->m_velo[21+i];
    v1467 = elemDat->m_velo[0+i] - elemDat->m_velo[9+i] - elemDat->m_velo[15+i] + elemDat->m_velo[18+i];
    v1256 = elemDat->m_velo[0+i] - elemDat->m_velo[3+i] - elemDat->m_velo[12+i] + elemDat->m_velo[15+i];

    hgr[i][0] = v1467 - v2358;
    hgr[i][1] = v1467 + v2358;
    hgr[i][2] = v1256 - v3478;
    hgr[i][3] = v1256 + v3478;
  }

  double ah = c, al = 100. * para;
  for(int i = 0; i < 3; i++){
    for(int k = 0; k < 4; k++){
      hgr[i][k] = hgr[i][k] * (ah + abs(al*hgr[i][k]));
    }
  }
  double hap, ham, hbp, hbm;
  temp.resize(elemDat->m_fint.size(), 0.);
  for(int i = 0; i < 3; i++){
    hap = hgr[i][0] + hgr[i][3];
    ham = hgr[i][0] - hgr[i][3];
    hbp = hgr[i][1] + hgr[i][2];
    hbm = hgr[i][1] - hgr[i][2];

    temp[ 0+i] = -hap - hbp;
    temp[ 3+i] =  hap - hbm;
    temp[ 6+i] = -hap + hbp;
		temp[ 9+i] =  hap + hbm;
		temp[12+i] = -ham + hbp;
		temp[15+i] =  ham + hbm;
		temp[18+i] = -ham - hbp;
		temp[21+i] =  ham - hbm;
  }

  return temp;
}