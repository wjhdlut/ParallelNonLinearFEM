
#ifndef OUTPUTWRITER_H
#define OUTPUTWRITER_H

#include <util/BaseModule.h>
#include <util/ObjectFactory.h>

class OutputWriter : public BaseModule
{
public:
  /**
   * @Brief: Construct a new Output Writer object
   * 
   * @param props 
   */
  OutputWriter(const nlohmann::json&props);

  /**
   * @Brief: Destroy the Output Writer object
   * 
   */
  ~OutputWriter();

  /**
   * @Brief:
   * 
   */
  virtual void Run() override;

private:
  bool m_onScreen = false;
};

ReflectRegister(OutputWriter, const nlohmann::json&)

#endif // OUTPUTWRITER_H