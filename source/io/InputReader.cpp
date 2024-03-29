
#include <io/InputReader.h>
#include <util/FileParse.h>
#include <fem/NodeSet.h>
#include <fem/ElementSet.h>
#include <fem/DofSpace.h>

namespace NONLINEARFEMIO
{
GlobalData* InputReader(int rank, char **args)
{
  std::string fileName = args[1];
  
  nlohmann::json props = nlohmann::json::object();
  FileParse(props, fileName);

  std::string dataFileName = props.at("input");
  std::shared_ptr<NodeSet> nodes = std::make_shared<NodeSet>();
  nodes->ReadFromFile(dataFileName);

  std::shared_ptr<ElementSet> elems = std::make_shared<ElementSet>(nodes, props);
  elems->ReadFromFile(dataFileName);
  
  std::shared_ptr<DofSpace> dofs = std::make_shared<DofSpace>(elems, nodes);
  dofs->ReadFromFile(dataFileName);

  GlobalData *globalData = GlobalData::GetInstance();
  globalData->SetFEMData(props, nodes, elems, dofs);
  globalData->ReadFromFile(dataFileName);

  GlobalData::GetInstance()->m_prefix = fileName.substr(0, fileName.find_first_of("."));

  return globalData;
}
}