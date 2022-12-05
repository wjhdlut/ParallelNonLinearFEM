#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <vector>

#include <util/BaseModule.h>

class OutputManager
{
public:
  OutputManager();
  ~OutputManager();

  void Run();

private:
  std::vector<std::shared_ptr<BaseModule>> m_outman;
};

#endif // OUTPUTMANAGER_H