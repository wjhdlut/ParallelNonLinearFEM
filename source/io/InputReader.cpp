
#include <io/InputReader.h>
#include <util/FileParse.h>
#include <fem/NodeSet.h>

namespace NONLINEARFEMIO
{
void InputReader(int rank)
{
  std::string fileName = "/home/wangjianhua/ExampleCpp/ParallelNonLinearFEM/input/cantilever8.pro";
  nlohmann::json bd = nlohmann::json::object();
  FileParse(bd, fileName);

  std::string dataFileName = bd.at("input");
  std::shared_ptr<NodeSet> nodes = std::make_shared<NodeSet>();
  nodes->ReadFromFile(dataFileName);
}
}