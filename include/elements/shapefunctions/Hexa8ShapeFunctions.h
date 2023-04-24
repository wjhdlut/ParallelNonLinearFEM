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
  virtual void GetShapeFunction(const std::vector<double> &xi) override;

  virtual std::vector<double> HourGlassTech(std::shared_ptr<ElementData> &elemDat,
                                            const double &c,
                                            const nlohmann::json &hourGlassPara,
                                            const std::vector<std::vector<double>> &pHpX) override;
  
private:
  std::vector<double> HourGlassFlangan(std::shared_ptr<ElementData> &elemDat,
                                       const double &c,
                                       const std::vector<std::vector<double>> &pHpX);
  
  std::vector<double> HourGlassStand(std::shared_ptr<ElementData> &elemDat,
                                     const double &c,
                                     const double &para);

private:
  std::vector<std::vector<double>> beat;
  std::vector<std::vector<double>> hgr;
  std::vector<std::vector<double>> m_h = {{1, -1, 1, -1},
                                          {1, 1, -1, -1},
                                          {1, -1, -1, 1},
                                          {1, -1, 1, -1}};
                                          
  std::vector<std::vector<double>> m_ss = {{2., -2., 2., -2.},
                                           {-2., -2., 2., 2.},
                                           {-2., 2., 2., -2.},
                                           {0., 0., 0., 0.}};
};

ReflectRegister(Hexa8ShapeFunctions)

#endif // HEXA8SHAPEFUNCTION_H