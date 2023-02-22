#ifndef TRUSS_H
#define TRUSS_H

#include <elements/Element.h>
#include <util/ObjectFactory.h>

class Truss : public Element
{
public:
  /**
   * @Brief: Construct a new Truss object
   * 
   * @param elemNode 
   * @param modelProps 
   */
  Truss(const std::vector<int> &elemNode, const nlohmann::json &modelProps);
  
  /**
   * @Brief: Destroy the Truss object
   * 
   */
  ~Truss();
};


#endif //TRUSS_H