
#include <util/ShapeFunctions.h>
#include <util/ObjectFactory.h>
#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/Math.h>
#include <iostream>

namespace ShapeFunctions
{

std::string GetElemType(const MatrixXd&elemCoords)
{
  int numOfNode = elemCoords.rows();
  int numOfDim  = elemCoords.cols();

  // 1-d element
  if(1 == numOfDim){
    if(2 == numOfNode){
      return "Line2";
    }
    else if(3 == numOfNode){
      return "Line3";
    }
    else{
      throw "No 1D element with " + std::to_string(numOfNode) + " nodes available";
    }
  }
  else if(2 == numOfDim){
    // 2-d element
    if(3 == numOfNode){
      return "Tria3";
    }
    else if(4 == numOfNode){
      return "Quad4";
    }
    else if(6 == numOfNode){
      return "Tria6";
    }
    else if(8 == numOfNode){
      return "Quad8";
    }
    else if(9 == numOfNode){
      return "Quad9";
    }
    else{
      throw "No 2D element with " + std::to_string(numOfNode) + " nodes available";
    }
  }
  else if(3 == numOfDim){
    // 3-d element
    if(4 == numOfNode){
      return "Tetra4";
    }
    else if(6 == numOfNode){
      return "Penta6";
    }
    else if(8 == numOfNode){
      return "Hexa8";
    }
    else{
      throw "No 3D element with " + std::to_string(numOfNode) + " nodes available";
    }
  }
  else{
    throw "the dimension must be 1,2 or 3";
  }
}

void GaussScheme(VectorXd&xi, VectorXd&weight, const int order)
{
  xi.resize(order), weight.resize(order);
  if(1 == order){
    xi(0) = 0., weight(0) = 2.;
  }
  if(2 == order){
    xi(0) = -0.577350269189626, weight(0) = 1.;
    xi(1) = +0.577350269189626, weight(1) = 1.;
  }
  if(3 == order){
    xi(0) = -0.774596669241483, weight(0) = 0.555555555555556;
    xi(1) = +0.,                weight(1) = 0.888888888888889;
    xi(2) = +0.774596669241483, weight(2) = 0.555555555555556;
  }
}

void TriaScheme(MatrixXd&xi, VectorXd&weight, const int order)
{
  xi.resize(order, 2), weight.resize(order);
  if(1 == order){
    xi(0, 0) = 1./3., xi(0, 1) = 1./3.;
    weight(0) = 1.0;
  }
  else if(3 == order){
    double r1 = 1./6., r2 = 2./3.;
    xi(0, 0) = r1, xi(0, 1) = r1;
    xi(1, 0) = r2, xi(1, 1) = r1;
    xi(2, 0) = r1, xi(2, 1) = r2;
    weight << 1./3., 1./3., 1./3.;
  }
  else if(7 == order){
    double r1 = 0.1012865073235;
    double r2 = 0.7974269853531;
    double r4 = 0.4701420641051;
    double r6 = 0.0597158717898;
    double r7 = 1.0/3.0;
    // xi = {{r1, r1}, {r2, r1}, {r1, r2}, {r4, r6}, {r4, r4}, {r6, r4}, {r7, r7}};
    xi(0, 0) = r1, xi(0, 1) = r1;
    xi(1, 0) = r2, xi(1, 1) = r1;
    xi(2, 0) = r1, xi(2, 1) = r2;
    xi(3, 0) = r4, xi(3, 1) = r6;
    xi(4, 0) = r4, xi(4, 1) = r4;
    xi(5, 0) = r6, xi(5, 1) = r4;
    xi(6, 0) = r7, xi(6, 1) = r7;

    double w1 = 0.1259391805448;
    double w4 = 0.1323941527885;
    double w7 = 0.225;
    weight << w1, w1, w1, w4, w4, w4, w7;
  }
  
}

void GetIntegrationPoints(MatrixXd&xi,
                          VectorXd&weight, 
                          const std::string&elemType,
                          const int order)
{
  if(0 != xi.size()) return;
  int stdOrder = 0;
  VectorXd xi1d, weight1d;
  if(elemType.npos != elemType.find("Line")){
    // 1-d element
    if("Line2" == elemType) stdOrder = 2;
    if("Line3" == elemType) stdOrder = 3;
    GaussScheme(xi1d, weight1d, stdOrder + order);
    xi = MatrixXd::Zero(xi1d.size(), 1);
    xi.col(0) = xi1d;
    weight = weight1d;
  }
  else if(elemType.npos != elemType.find("Tria")){
    // 2-d triangle element
    std::vector<int> orderArray = {1, 3, 7};
    if("Tria3" == elemType) stdOrder = 0;
    if("Tria3" == elemType) stdOrder = 1;
    TriaScheme(xi, weight, orderArray[stdOrder + order]);
  }
  else if(elemType.npos != elemType.find("Quad")){
    // 2-d quadrilateral element 
    if("Quad4" == elemType) stdOrder = 2;
    if( "Quad8" == elemType || "Quad9" == elemType) stdOrder = 3;
    GaussScheme(xi1d, weight1d, stdOrder + order);
    xi.resize(std::pow(stdOrder + order, 2), 2);
    weight.resize(std::pow(stdOrder + order, 2));
    int index = 0;
    for(int i = 0; i < stdOrder + order; i++){
      for(int j = 0; j < stdOrder + order; j++){
        xi(index, 0) = xi1d(i), xi(index, 1) = xi1d(j);
        weight(index) = (weight1d(i) * weight1d(j));
        index += 1;
      }
    }
  }
  else if(elemType.npos != elemType.find("Hexa")){
    // 3-d hexahedron element
    if("Hexa8" == elemType) stdOrder = 2;
    GaussScheme(xi1d, weight1d, stdOrder + order);
    xi.resize(std::pow(stdOrder + order, 3), 3);
    weight.resize(std::pow(stdOrder + order, 3));
    int index = 0;
    for(int i = 0; i < stdOrder + order; i++){
      for(int j = 0; j < stdOrder + order; j++){
        for(int k = 0; k < stdOrder + order; k++){
          xi(index, 0) = xi1d(i), xi(index, 1) = xi1d(j), xi(index, 2) = xi1d(k);
          weight(index) = weight1d(i) * weight1d(j) * weight1d(k);
          index += 1;
        }
      }
    }
  }
}

void GetElemShapeData(MatrixXd&elemCoords, const int order,
                      const std::string &method, const std::string&elemType)
{
  // std::string realElemType = ("Default" == elemType) ? GetElemType(elemCoords) : elemType;
  // std::vector<std::vector<double>> xi;
  // std::vector<double> weight;
  // GetIntegrationPoints(xi, weight, realElemType, order, method);
  // std::string elemName = realElemType + "ShapeFunctions";

  // std::vector<std::vector<double>> jac;
  // std::vector<std::vector<double>> invJac;
  
  // /* pHpX = [pH1pX1 pH1pX2 pH1pX3,
  //            pH2pX1 pH2pX2 pH2pX3,
  //            ...
  //            pHnpX1 pHnpX2 pHnpX3,]*/
  // std::vector<std::vector<double>> pHpX;

  // double weighti = 0.;
  // int count = 0;
  // for(auto iXi : xi){
  //   std::shared_ptr<ElementShapeFunctions>res = ObjectFactory::CreateObject<ElementShapeFunctions>(elemName);
  //   if(nullptr == res) throw "Unkonwn type " + realElemType;
  //   res->GetShapeFunction(iXi);

  //   // compute jacobian matrix
  //   /*jac = [pXpxi1 pXpxi2 pXpxi3,
  //            pYpxi1 pYpxi2 pYpxi3,
  //            pZpxi1 pZpxi2 pZpxi3];*/
  //   jac = Math::MatrixATransMultB(elemCoords, res->pHpxi);

  //   invJac = Math::MatrixInverse(jac);

  //   pHpX = Math::MatrixAMultB(res->pHpxi, invJac);

  //   weighti = Math::MatrixDet(jac) * weight[count];

  //   count++;
  // }
}
}