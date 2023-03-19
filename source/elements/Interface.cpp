#include <elements/Interface.h>
#include <util/ShapeFunctions.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>

Interface::Interface(const std::vector<int> &elemNodes, const nlohmann::json &modelProps)
          : Element(elemNodes, modelProps)
{
  method = "NewtonCotes";

  std::vector<double> temp(2, 0.);
  SetHistoryParameter("normal", temp);

  CommitHistory();
}

Interface::~Interface()
{}

void Interface::GetTangentStiffness(std::shared_ptr<ElementData> &elemDat)
{
  Matrix rot = GetRotation(elemDat->m_coords, elemDat->m_state);

  std::string elemType = ShapeFunctions::GetElemType(elemDat->m_coords);
  ShapeFunctions::GetIntegrationPoints(xi, weight, elemType, order, method);

  std::string elemName = elemType + "ShapeFunctions";
  std::shared_ptr<ElementShapeFunctions> res = ObjectFactory::CreateObject<ElementShapeFunctions>(elemName);
  if(nullptr == res) throw "Unknown type " + elemType;

  elemDat->m_outLabel.emplace_back("tractions");
  outputData = Math::MatrixZeros(elemDat->m_coords.size(), 2);

  int count = 0;
  Matrix tempMatrix;
  for(auto iXi : xi)
  {
     res->GetShapeFunction(iXi);

     GetBMatrix(res->H, rot);

     GetKinematics(elemDat->m_state);

     sigma = m_mat->GetStress(kin);

     D = m_mat->GetTangMatrix();

     tempMatrix = Math::MatrixAMultB(D, B);
     elemDat->m_stiff = Math::MatrixATransMultB(B, tempMatrix);
     elemDat->m_stiff = Math::MatrixScale(weight[count], elemDat->m_stiff);

     elemDat->m_fint = Math::MatrixATransMultVecB(B, sigma);
     elemDat->m_fint = Math::VecScale(weight[count], elemDat->m_fint);

     tempMatrix = Math::VecOuter(std::vector<double>(elemDat->m_coords.size(), 1.), sigma);
     elemDat->m_outputData = Math::MatrixAdd(1., elemDat->m_outputData, tempMatrix);
  }
}

Matrix Interface::GetRotation(const Matrix &coords, const std::vector<double> &state)
{
  Matrix midCoords = Math::MatrixZeros(2, 2);
  
  for(int i = 0; i < 2; i++)
    midCoords[i] = Math::VecAdd(0.5, coords[i], coords[2-i]);

  midCoords[0][0] += 0.5 * ( state[0] + state[4] );
  midCoords[0][1] += 0.5 * ( state[1] + state[5] );
  midCoords[1][0] += 0.5 * ( state[2] + state[6] );
  midCoords[1][1] += 0.5 * ( state[3] + state[7] );

  std::vector<double> ds = Math::VecAdd(-1., midCoords[1], midCoords[0]);

  std::vector<double> normal = GetHistoryParameter("normal");
  
  if(Math::VecNorm(normal) < 0.5){
    normal[0] = ds[1] / Math::VecNorm(ds);
    normal[1] = ds[0] / Math::VecNorm(ds);
  }
  else{
    std::vector<double> newNormal(normal.size(), 0.);
    newNormal[0] = ds[1] / Math::VecNorm(ds);
    newNormal[1] = ds[0] / Math::VecNorm(ds);

    if(Math::VecDot(newNormal, normal) < 0.)
      Math::VecScale(-1., newNormal);
    
    normal.swap(newNormal);
  }
  SetHistoryParameter("normal", normal);

  Matrix rot = {{normal[0], normal[1]}, {normal[1], normal[0]}};

  return rot;
}

void Interface::GetBMatrix(const std::vector<double> &H, const Matrix &R)
{
  int numOfDim = m_dofType.size();
  int numOfNode = H.size();
  std::vector<double> temp(numOfDim*numOfNode, 0.);
  B.resize(2, temp);

  B[0][0] = -R[0][0] * H[0], B[0][1] = -R[0][1] * H[0];
  B[1][0] = -R[1][0] * H[0], B[1][1] = -R[1][1] * H[0];

  B[0][2] = -R[0][0] * H[1], B[0][3] = -R[0][1] * H[1];
  B[1][2] = -R[1][0] * H[1], B[1][3] = -R[1][1] * H[1];

  B[0][4] = R[0][0] * H[0], B[0][5] = R[0][1] * H[0];
  B[1][4] = R[1][0] * H[0], B[1][5] = R[1][1] * H[0];

  B[0][6] = R[0][0] * H[1], B[0][7] = R[0][1] * H[1];
  B[1][6] = R[1][0] * H[1], B[1][7] = R[1][1] * H[1];
}

void Interface::GetKinematics(const std::vector<double> &elState)
{
  int numOfDim = 2;
  
  kin = std::make_shared<Kinematics>(numOfDim);

  kin->strain = Math::MatrixAMultVecB(B, elState);
}