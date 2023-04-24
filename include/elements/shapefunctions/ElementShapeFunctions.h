
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

  virtual std::vector<double> HourGlassTech(std::shared_ptr<ElementData> &elemDat,
                                            const double &c,
                                            const nlohmann::json &para,
                                            const std::vector<std::vector<double>> &pHpX)
  {
    return std::vector<double>(0.);
  }
  
  inline std::vector<std::vector<int>> *ReturnFaceMatrix(){
    return &m_face;
  }

public:
  int numOfStress;                              // The number of stress
  std::vector<double> H;                        // The shape function
  std::vector<std::vector<double>> pHpxi;       // The derivative of shape function about local coordinate system
  std::vector<std::string> dofType;             // Element Dof Type

protected:
  std::vector<std::vector<int>> m_face;
};

#endif // ELEMENTSHAPEFUNCTIONS_H
