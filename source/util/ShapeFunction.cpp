
#include <util/ShapeFunctions.h>
#include <util/ObjectFactory.h>
#include <elements/ElementShapeFunctions.h>
#include <util/Math.h>

namespace ShapeFunctions
{

std::string GetElemType(const std::vector<std::vector<double>>&elemCoords)
{
  int numOfNode = elemCoords.size();
  int numOfDim = elemCoords[0].size();

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

void GaussScheme(std::vector<double>&xi, std::vector<double>&weight, const int order)
{
  if(2 == order){
    xi.emplace_back(-0.577350269189626), weight.emplace_back(1.);
    xi.emplace_back(+0.577350269189626), weight.emplace_back(1.);
  }
  if(3 == order){
    xi.emplace_back(-0.774596669241483), weight.emplace_back(0.555555555555556);
    xi.emplace_back(+0.),                weight.emplace_back(0.888888888888889);
    xi.emplace_back(+0.774596669241483), weight.emplace_back(0.555555555555556);
  }
}

void TriaScheme(std::vector<std::vector<double>>&xi, std::vector<double>&weight, const int order)
{
  if(1 == order){
    xi = {{1./3., 1./3.}};
    weight.emplace_back(1.0);
  }
  else if(3 == order){
    double r1 = 1./6., r2 = 2./3.;
    xi = {{r1, r1}, {r2, r1}, {r1, r2}};
    weight = {1./3., 1./3., 1./3.};
  }
  else if(7 == order){
    double r1 = 0.1012865073235;
    double r2 = 0.7974269853531;
    double r4 = 0.4701420641051;
    double r6 = 0.0597158717898;
    double r7 = 1.0/3.0;
    xi = {{r1, r1}, {r2, r1}, {r1, r2}, {r4, r6}, {r4, r4}, {r6, r4}, {r7, r7}};

    double w1 = 0.1259391805448;
    double w4 = 0.1323941527885;
    double w7 = 0.225;
    weight = {w1, w1, w1, w4, w4, w4, w7};
  }
  
}

void GetIntegrationPoints(std::vector<std::vector<double>>&xi, std::vector<double>&weight,
                          const std::string&elemType, const int order, const::std::string &method)
{
  int stdOrder = 0;
  std::vector<double> xi1d;
  std::vector<double> weight1d;
  if(elemType.npos != elemType.find("Line")){
    // 1-d element
    if("Line2" == elemType) stdOrder = 2;
    if("Line3" == elemType) stdOrder = 3;
    GaussScheme(xi1d, weight1d, stdOrder + order);
    xi.emplace_back(xi1d);
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
    for(int i = 0; i < stdOrder + order; i++){
      for(int j = 0; j < stdOrder + order; j++){
        std::vector<double> iXi;
        iXi.emplace_back(xi1d[i]);
        iXi.emplace_back(xi1d[j]);
        xi.emplace_back(iXi);
        weight.emplace_back(weight1d[i] * weight1d[j]);
      }
    }
  }
  else if(elemType.npos != elemType.find("Hexa")){
    // 3-d hexahedron element
    if("Hexa8" == elemType) stdOrder = 2;
    GaussScheme(xi1d, weight1d, stdOrder + order);
    for(int i = 0; i < stdOrder + order; i++){
      for(int j = 0; j < stdOrder + order; j++){
        for(int k = 0; k < stdOrder + order; k++){
          std::vector<double> iXi;
          iXi.emplace_back(xi1d[i]);
          iXi.emplace_back(xi1d[j]);
          iXi.emplace_back(xi1d[k]);
          xi.emplace_back(iXi);
          weight.emplace_back(weight1d[i] * weight1d[j] * weight1d[k]);
        }
      }
    }
  }
}

void GetElemShapeData(std::vector<std::vector<double>>&elemCoords, const int order,
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