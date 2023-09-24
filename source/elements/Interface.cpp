#include <elements/Interface.h>
#include <util/ShapeFunctions.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>

Interface::Interface(const std::vector<int> &elemNodes, const nlohmann::json &modelProps)
          : Element(elemNodes, modelProps)
{
  method = "NewtonCotes";

  Vector2d temp = Vector2d::Zero();
  SetHistoryParameter("normal", temp);

  CommitHistory();
}

Interface::~Interface()
{}

void Interface::GetTangentStiffness(std::shared_ptr<ElementData> &elemDat)
{
  MatrixXd rot = GetRotation(elemDat->m_coords, elemDat->m_state);

  std::string elemType = ShapeFunctions::GetElemType(elemDat->m_coords);
  ShapeFunctions::GetIntegrationPoints(xi, weight, elemType, order, method);

  std::string elemName = elemType + "ShapeFunctions";
  std::shared_ptr<ElementShapeFunctions> res
              = ObjectFactory::CreateObject<ElementShapeFunctions>(elemName);
  if(nullptr == res) throw "Unknown type " + elemType;

  elemDat->m_outLabel.emplace_back("tractions");
  outputData = MatrixXd::Zero(elemDat->m_coords.rows(), 2);

  int count = 0;
  MatrixXd tempMatrix;
  for(int i = 0; i < xi.rows(); i++)
  {
    res->GetShapeFunction(xi.row(i));

    GetBMatrix(res->H, rot);

    GetKinematics(elemDat->m_state);

    sigma = m_mat->GetStress(kin, elemDat->m_Dstate);

    D = m_mat->GetTangMatrix();
    
    elemDat->m_stiff += weight[i] * (B.transpose() * D * B);

    elemDat->m_fint += weight[i] * (B.transpose() * sigma);
    
    elemDat->m_outputData += Math::VecCross(VectorXd::Ones(elemDat->m_coords.rows()), sigma);
  }
}

MatrixXd Interface::GetRotation(const MatrixXd &coords, const VectorXd &state)
{
  Matrix2d midCoords = Matrix2d::Zero();
  
  for(int i = 0; i < 2; i++)
    midCoords.row(i) = 0.5 * (coords.row(i) + coords.row(2-i));

  midCoords(0, 0) += 0.5 * ( state(0) + state(4) );
  midCoords(0, 1) += 0.5 * ( state(1) + state(5) );
  midCoords(1, 0) += 0.5 * ( state(2) + state(6) );
  midCoords(1, 1) += 0.5 * ( state(3) + state(7) );

  Vector2d ds = midCoords.row(1) - midCoords.row(0);

  Vector2d normal = GetHistoryParameter("normal");
  
  if(normal.norm() < 0.5){
    normal(0) = ds(1) / ds.norm();
    normal(1) = ds(0) / ds.norm();
  }
  else{
    Vector2d newNormal(0., 0.);
    newNormal[0] = ds[1] / ds.norm();
    newNormal[1] = ds[0] / ds.norm();

    if(newNormal.dot(normal) < 0.)
      newNormal = -1. * newNormal;
    
    normal = newNormal;
  }
  SetHistoryParameter("normal", normal);

  Matrix2d rot;
  rot << normal[0], normal[1], normal[1], normal[0];

  return rot;
}

void Interface::GetBMatrix(const VectorXd &H, const MatrixXd &R)
{
  int numOfDim = m_dofType.size();
  int numOfNode = H.rows();
  B = MatrixXd::Zero(2, numOfDim*numOfNode);

  B(0, 0) = -R(0, 0) * H(0), B(0, 1) = -R(0, 1) * H(0);
  B(1, 0) = -R(1, 0) * H(0), B(1, 1) = -R(1, 1) * H(0);

  B(0, 2) = -R(0, 0) * H(1), B(0, 3) = -R(0, 1) * H(1);
  B(1, 2) = -R(1, 0) * H(1), B(1, 3) = -R(1, 1) * H(1);

  B(0, 4) = R(0, 0) * H(0), B(0, 5) = R(0, 1) * H(0);
  B(1, 4) = R(1, 0) * H(0), B(1, 5) = R(1, 1) * H(0);

  B(0, 6) = R(0, 0) * H(1), B(0, 7) = R(0, 1) * H(1);
  B(1, 6) = R(1, 0) * H(1), B(1, 7) = R(1, 1) * H(1);
}

void Interface::GetKinematics(const VectorXd &elState)
{
  int numOfDim = 2;
  kin = std::make_shared<Kinematics>(numOfDim);
  kin->strain = B.transpose() * elState;
}