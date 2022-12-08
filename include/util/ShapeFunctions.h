#ifndef SHAPEFUNCTIONS_H
#define SHAPEFUNCTIONS_H

#include <vector>
#include <string>

namespace ShapeFunctions
{
void GaussScheme(std::vector<double>&xi, std::vector<double>&weight, const int order);

void TriaScheme(const int order);

std::string GetElemType(const std::vector<std::vector<double>>&elemCoords);

void GetIntegrationPoints(std::vector<std::vector<double>>&xi, std::vector<double>&weight,
                          const std::string&elemType, const int order, const::std::string &method);

void GetElemShapeData(std::vector<std::vector<double>>&elemCoords, const int order = 0,
                      const std::string &method = "Gauss", const std::string&elemType = "Default");
}

#endif // SHAPEFUNCTIONS_H
