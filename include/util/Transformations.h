#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <vector>


namespace Transformations
{
std::vector<double> ToElementCoordinates(const std::vector<double> &a,
                                         const std::vector<std::vector<double>> &elemNodeCoords);

std::vector<std::vector<double>> ToElementCoordinates(const std::vector<std::vector<double>> &A,
                                                      const std::vector<std::vector<double>> &elemNodeCoords);

std::vector<double> ToGlobalCoordinates(const std::vector<double> &a,
                                        const std::vector<std::vector<double>> &elemNodeCoords);

std::vector<std::vector<double>> ToGlobalCoordinates(const std::vector<std::vector<double>> &A,
                                                     const std::vector<std::vector<double>> &elemNodeCoords);

std::vector<double> VecToElementCoordinates(const std::vector<double> &a,
                                            const std::vector<std::vector<double>> &elemNodeCoords);

std::vector<std::vector<double>> MatToElementCoordinates(const std::vector<std::vector<double>> &A,
                                                         const std::vector<std::vector<double>> &elemNodeCoords);

std::vector<double> VecToGLobalCoordinates(const std::vector<double> &a,
                                           const std::vector<std::vector<double>> &elemNodeCoords);

std::vector<std::vector<double>> MatToGlobalCoordinates(const std::vector<std::vector<double>> &A,
                                                        const std::vector<std::vector<double>> &elemNodeCoords);

std::vector<std::vector<double>> GetRotationMatrix(const std::vector<std::vector<double>> &elemNodeCoords);

}

#endif // TRANSFORMATIONS_H