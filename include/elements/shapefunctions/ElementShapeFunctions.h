
#ifndef ELEMENTSHAPEFUNCTIONS_H
#define ELEMENTSHAPEFUNCTIONS_H

#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include <elements/Element.h>

/**
 * @Brief:   Base Class to Compute Element Shape Functions
 * 
 */

/* pHpxi = [pH1pxi1 pH1pxi2 pH1pxi3,
            pH2pxi1 pH2pxi2 pH2pxi3，
            ...,
            pHnpxi1 pHnpxi2 pHnpxi3]*/

struct ElementData;

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
  virtual void GetShapeFunction(const VectorXd &xi) = 0;

  virtual VectorXd HourGlassTech(std::shared_ptr<ElementData> &elemDat,
                                 const VectorXd &elemNodeDisp,
                                 const double &c,
                                 const nlohmann::json &para,
                                 const MatrixXd &pHpX)
  {
    return VectorXd::Zero(0);
  }
  
  inline MatrixXi *ReturnFaceMatrix(){
    return &m_face;
  }

  virtual double ComputeElemTimeStep(double &dtK1, double &elemDistortion,
                                     const MatrixXd &elemNodeCoords,
                                     const VectorXd &elemNodeDisp,
                                     const double detJac, const double waveSpeed)
  {
    return 0.;
  }

public:
  int numOfStress;                              // The number of stress
  VectorXd H;                                   // The shape function
  MatrixXd pHpxi;                               // The derivative of shape function about local coordinate system
  std::vector<std::string> dofType;             // Element Dof Type

protected:
  MatrixXi m_face;
  std::vector<std::string> m_dofType;
};

#endif // ELEMENTSHAPEFUNCTIONS_H
