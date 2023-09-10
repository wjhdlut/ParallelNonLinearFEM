#ifndef HEXA8SHAPEFUNCTION_H
#define HEXA8SHAPEFUNCTION_H

#include <elements/shapefunctions/ElementShapeFunctions.h>
#include <util/ObjectFactory.h>

/* ----------------------------------------------------
 *  element shape and node order
 *                    5 -------------------- 8
 *                   / |                    /|
 *                  /  |                   / |
 *                 /   |                  /  |
 *                /    |                 /   |
 *               /     |                /    |
 *              6 -------------------- 7     |
 *              |      |               |     |
 *              |      1 --------------|---- 4
 *              |     /                |     /
 *              |    /                 |    /
 *              |   /                  |   /
 *              |  /                   |  /
 *              | /                    | /
 *              2 ---------------------3
 * ----------------------------------------------------- */ 

class Hexa8ShapeFunctions : public ElementShapeFunctions
{
public:
  /**
   * @Brief: Construct Shape Functions object for 8 Nodes Hexahedron
   * 
   */
  Hexa8ShapeFunctions();

  /**
   * @Brief: Destroy Shape Functions object for 8 Nodes Hexahedron
   * 
   */
  ~Hexa8ShapeFunctions();

  /**
   * @Brief:  Compute Shape Functions for 8 Nodes Hexahedron
   * 
   * @param xi     [in]  Gauss Points Coordinates
   */
  virtual void GetShapeFunction(const VectorXd &xi) override;

  virtual VectorXd HourGlassTech(std::shared_ptr<ElementData> &elemDat,
                                 const VectorXd &elemNodeDisp, 
                                 const double &c,
                                 const nlohmann::json &hourGlassPara,
                                 const MatrixXd &pHpX) override;
  
  virtual double ComputeElemTimeStep(double &dtK1, double &elemDistortion,
                                   const MatrixXd &elemNodeCoords,
                                   const VectorXd &elemNodeDisp,
                                   const double detJac, const double waveSpeed) override;
private:
  VectorXd HourGlassFlangan(std::shared_ptr<ElementData> &elemDat,
                            const VectorXd &elemNodeDisp,
                            const double &c,
                            const MatrixXd &pHpX);
  
  VectorXd HourGlassStand(std::shared_ptr<ElementData> &elemDat,
                          const double &c,
                          const double &para);

private:
  MatrixXd beat;
  MatrixXd hgr;
  Matrix4d m_h;
                                          
  Matrix4d m_ss;
};

ReflectRegister(Hexa8ShapeFunctions)

#endif // HEXA8SHAPEFUNCTION_H