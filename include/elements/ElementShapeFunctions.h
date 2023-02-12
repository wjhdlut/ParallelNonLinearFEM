
#ifndef ELEMENTSHAPEFUNCTIONS_H
#define ELEMENTSHAPEFUNCTIONS_H

#include <vector>

/**
 * @Brief:   Base Class to Compute Element Shape Functions
 * 
 */

/* pHpxi = [pH1pxi1 pH1pxi2 pH1pxi3,
            pH2pxi1 pH2pxi2 pH2pxi3，
            ...,
            pHnpxi1 pHnpxi2 pHnpxi3]*/

class ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct a new Element Shape Functions object
   * 
   */
  ElementShapeFunctions(){}

  /**
   * @Brief: Destroy the Element Shape Functions object
   * 
   */
  virtual ~ElementShapeFunctions(){}

  /**
   * @Brief: Compute the Shape Function Virtual Method
   * 
   * @param xi                [in]  Gauss Coordinate
   */
  virtual void GetShapeFunction(const std::vector<double> &xi) = 0;

public:
  int numOfStress;                              // The number of stress
  std::vector<double> H;                        // the shape function
  std::vector<std::vector<double>> pHpxi;       // the derivative of shape function about local coordinate system
};

#endif // ELEMENTSHAPEFUNCTIONS_H
