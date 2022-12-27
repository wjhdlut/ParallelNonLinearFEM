#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <vector>
#include <util/Math.h>

class Kinematics{
public:
  inline Kinematics(const int dim){
    F = Math::MatrixEye(dim);
    
    std::vector<double> rowData(dim, 0);
    E.resize(dim, rowData);

    int numOfStrain = dim * (dim + 1)/2;
    strain.resize(numOfStrain, 0);

    numOfDim = dim;
  }

  ~Kinematics(){}

  inline void SetStrainVector(){
    if(2 == numOfDim){
      strain[0] = E[0][0];
      strain[1] = E[1][1];
      strain[2] = E[0][1];
    }
    else if(3 == numOfDim){
      strain[0] = E[0][0];
      strain[1] = E[1][1];
      strain[2] = E[2][2];

      strain[3] = 2. * E[0][1];
      strain[4] = 2. * E[1][2];
      strain[5] = 2. * E[0][2];
    }
  }

public:
  std::vector<std::vector<double>> F;
  std::vector<std::vector<double>> E;
  std::vector<double> strain;

private:
  int numOfDim = 0;
};

#endif // KINEMATICS_H

