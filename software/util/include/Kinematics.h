/**
 * @File Name:     Kinematics.h
 * @Author:        JianHuaWang (992411152@qq.com)
 * @Brief:         
 * @Version:       0.1
 * @Create Date:   2023-10-25
 * 
 * @Copyright Copyright (c) 2023 JianHuaWang
 * 
 */

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
    incremStrain = VectorXd::Zero(numOfStrain);

    numOfDim = dim;
  }

  ~Kinematics(){}

  /**
   * @Brief:   Set the Strain Vector object based on the Green Strain Tensor
   * 
   */
  inline void SetStrainVector(const MatrixXd &tensor){
    if(2 == numOfDim){
      strain(0) = tensor(0, 0);
      strain(1) = tensor(1, 1);
      strain(2) = 2. * tensor(0, 1);
    }
    else if(3 == numOfDim){
      strain(0) = tensor(0, 0);
      strain(1) = tensor(1, 1);
      strain(2) = tensor(2, 2);

      strain(3) = 2. * tensor(0, 1);
      strain(4) = 2. * tensor(1, 2);
      strain(5) = 2. * tensor(0, 2);
    }
  }

public:
  MatrixXd F;                   // Deformation Gradient Tensor
  MatrixXd E;                   // Green Strain Tensor
  VectorXd strain;              // Strain Vector
  VectorXd incremStrain;        // Incremental Strain Vector

private:
  int numOfDim = 0;
};

#endif // KINEMATICS_H

