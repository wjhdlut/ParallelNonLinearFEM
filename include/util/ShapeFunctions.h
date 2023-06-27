#ifndef SHAPEFUNCTIONS_H
#define SHAPEFUNCTIONS_H

#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace ShapeFunctions
{
void GaussScheme(VectorXd&xi, VectorXd&weight, const int order);

void TriaScheme(const int order);

std::string GetElemType(const MatrixXd&elemCoords);

void GetIntegrationPoints(MatrixXd&xi, VectorXd&weight, const std::string&elemType,
                          const int order, const::std::string &method);

void GetElemShapeData(MatrixXd&elemCoords, const int order = 0,
                      const std::string &method = "Gauss",
                      const std::string&elemType = "Default");
}

#endif // SHAPEFUNCTIONS_H
