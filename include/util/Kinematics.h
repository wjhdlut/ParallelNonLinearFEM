#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <vector>
// #include <util/Math.h>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

/**
 * @Brief:  compute Deformation Gradient tensor F, Green Strain E tensor,
 *  and strain vector
 * 
 */

/** ------------------------------------
 *     F = [px/pX, px/pY, px/pZ 
 *          py/pX, py/pY, py/pZ,
 *          pz/pX, pz/pY, pz/pZ]
 * ------------------------------------- */
class Kinematics{
public:
  // Initial the Deformation Gradient and
  // Green Strain Tensor
  inline Kinematics(const int dim){
    F.setIdentity(dim, dim);
    E = MatrixXd::Zero(dim, dim);

    int numOfStrain = dim * (dim + 1)/2;
    strain = VectorXd::Zero(numOfStrain);

    numOfDim = dim;
  }

  ~Kinematics(){}

  /**
   * @Brief:   Set the Strain Vector object based on the Green Strain Tensor
   * 
   */
  inline void SetStrainVector(){
    if(2 == numOfDim){
      strain(0) = E(0, 0);
      strain(1) = E(1, 1);
      strain(2) = 2. * E(0, 1);
    }
    else if(3 == numOfDim){
      strain(0) = E(0, 0);
      strain(1) = E(1, 1);
      strain(2) = E(2, 2);

      strain(3) = 2. * E(0, 1);
      strain(4) = 2. * E(1, 2);
      strain(5) = 2. * E(0, 2);
    }
  }

public:
  MatrixXd F;                   // Deformation Gradient Tensor
  MatrixXd E;                   // Green Strain Tensor
  VectorXd strain;              // strain vector

private:
  int numOfDim = 0;
};

#endif // KINEMATICS_H

