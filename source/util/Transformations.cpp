#include <util/Transformations.h>
#include <util/Math.h>

namespace Transformations
{

std::vector<double> ToElementCoordinates(const std::vector<double> &a,
                                         const std::vector<std::vector<double>> &elemNodeCoords)
{
  return VecToElementCoordinates(a, elemNodeCoords);
}

std::vector<std::vector<double>> ToElementCoordinates(const std::vector<std::vector<double>> &A,
                                                      const std::vector<std::vector<double>> &elemNodeCoords)
{
  return MatToElementCoordinates(A, elemNodeCoords);
}

std::vector<double> ToGlobalCoordinates(const std::vector<double> &a,
                                        const std::vector<std::vector<double>> &elemNodeCoords)
{
  return VecToGLobalCoordinates(a, elemNodeCoords);
}

std::vector<std::vector<double>> ToGlobalCoordinates(const std::vector<std::vector<double>> &A,
                                                     const std::vector<std::vector<double>> &elemNodeCoords)
{
  return MatToGlobalCoordinates(A, elemNodeCoords);
}

std::vector<double> VecToElementCoordinates(const std::vector<double> &a,
                                            const std::vector<std::vector<double>> &elemNodeCoords)
{
  Matrix R = GetRotationMatrix(elemNodeCoords);

  std::vector<double> aBar;

  if(0 != a.size()%R.size()) throw "Vector does not have the right shape to be rotated";

  std::vector<double> tempVec(R.size(), 0.);
  std::vector<double> tempVecA;
  for(int i = 0; i < a.size()/R.size(); i++)
  {
    tempVecA.clear();
    tempVecA.insert(tempVecA.begin(), a.begin()+R.size()*i, a.begin()+R.size()*(i+1));
    tempVec = Math::MatrixAMultVecB(R, tempVecA);
    aBar.insert(aBar.end(), tempVec.begin(), tempVec.end());
  }

  return aBar;
}

std::vector<std::vector<double>> MatToElementCoordinates(const std::vector<std::vector<double>> &A,
                                                         const std::vector<std::vector<double>> &elemNodeCoords)
{
  Matrix R = GetRotationMatrix(elemNodeCoords);

  std::vector<double> tempVec(A[0].size(), 0.);
  Matrix aBar(A.size(), tempVec);

  Matrix tempMatA, tempMat;
  
  tempVec.clear();
  if(0 != A[0].size() % R.size() || 0 != A.size() % R.size())
    throw "Matrix does not have the right shape to be rotated";
  
  for(int i = 0; i < A[0].size()/R.size(); i++)
  {
    for(int j = 0; j < A.size()/R.size(); j++)
    {
      tempMat.clear();
      for(int jj = 0; jj < A.size(); jj++)
      {
        tempVec.insert(tempVec.begin(), A[j*R.size()+jj].begin()+R.size()*i, A[j*R.size()+jj].begin()+R.size()*(i+1));
        tempMatA.emplace_back(tempVec);
      }
      tempMat = Math::MatrixAMultB(R, tempMatA);
      tempMat = Math::MatrixAMultBTrans(tempMat, R);

      for(int ii = 0; ii < R.size(); ii++)
        for(int jj = 0; jj < R.size(); jj++)
          aBar[j*R.size() + ii][i*R.size() + jj] = tempMat[ii][jj];
    }
  }

  return aBar;
}

std::vector<std::vector<double>> GetRotationMatrix(const std::vector<std::vector<double>> &elemNodeCoords)
{
  if(0 == elemNodeCoords.size()) return std::vector<std::vector<double>>(0);

  if(2 != elemNodeCoords[0].size()) throw "Rotation matrix only implemented for 2D situation";


  // Compute the undeformed element length
  double elemLength = Math::VecNorm(Math::VecAdd(-1., elemNodeCoords[1], elemNodeCoords[0]));

  // rotate a globdal coordinate to an element coordinate
  double sin_alpha = (elemNodeCoords[1][1] - elemNodeCoords[0][1])/elemLength;
  double cos_alpha = (elemNodeCoords[1][0] - elemNodeCoords[0][0])/elemLength;
  
  Matrix A = {{cos_alpha, sin_alpha},{-sin_alpha, cos_alpha}};
  return A;
}

std::vector<double> VecToGLobalCoordinates(const std::vector<double> &a,
                                           const std::vector<std::vector<double>> &elemNodeCoords)
{
  Matrix R = GetRotationMatrix(elemNodeCoords);

  std::vector<double> aBar;

  if(0 != a.size()%R.size()) throw "Vector does not have the right shape to be rotated";

  std::vector<double> tempVec(R.size(), 0.);
  std::vector<double> tempVecA;

  for(int i = 0; i < a.size()/R.size(); i++)
  {
    tempVecA.clear();
    tempVecA.insert(tempVecA.begin(), a.begin()+R.size()*i, a.begin()+R.size()*(i+1));
    tempVec = Math::MatrixATransMultVecB(R, tempVecA);
    aBar.insert(aBar.end(), tempVec.begin(), tempVec.end());
  }

  return aBar;
}

std::vector<std::vector<double>> MatToGlobalCoordinates(const std::vector<std::vector<double>> &A,
                                                        const std::vector<std::vector<double>> &elemNodeCoords)
{
  Matrix R = GetRotationMatrix(elemNodeCoords);

  std::vector<double> tempVec(A[0].size(), 0.);
  std::vector<std::vector<double>> aBar(A.size(), tempVec);

  Matrix tempMatA, tempMat;
  
  tempVec.clear();
  if(0 != A[0].size() % R.size() || 0 != A.size() % R.size())
    throw "Matrix does not have the right shape to be rotated";
  
  for(int i = 0; i < A[0].size()/R.size(); i++)
  {
    for(int j = 0; j < A.size()/R.size(); j++)
    {
      tempMatA.clear();
      tempVec.clear();
      for(int jj = 0; jj < R.size(); jj++)
      {
        tempVec.insert(tempVec.begin(), A[j*R.size()+jj].begin()+R.size()*i, A[j*R.size()+jj].begin()+R.size()*(i+1));
        tempMatA.emplace_back(tempVec);
      }
      tempMat = Math::MatrixATransMultB(R, tempMatA);
      tempMat = Math::MatrixAMultB(tempMat, R);

      for(int ii = 0; ii < R.size(); ii++)
        for(int jj = 0; jj < R.size(); jj++)
          aBar[j*R.size() + ii][i*R.size() + jj] = tempMat[ii][jj];
    }
  }

  return aBar;
}
}