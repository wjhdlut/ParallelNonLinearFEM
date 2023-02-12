#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <vector>

#include <util/BaseModule.h>

class OutputManager
{
public:
  /**
   * @Brief: Construct a new Output Manager object
   * 
   */
  OutputManager();

  /**
   * @Brief: Destroy the Output Manager object
   * 
   */
  ~OutputManager();

  /**
   * @Brief: output related variables
   * 
   */
  void Run();

private:
  std::vector<std::shared_ptr<BaseModule>> m_outman;
};

#endif // OUTPUTMANAGER_H