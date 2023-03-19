
#include <io/OutputWriter.h>
#include <util/DataStructure.h>

OutputWriter::OutputWriter(const nlohmann::json &props) : BaseModule(props)
{
  
}

OutputWriter::~OutputWriter()
{}

void OutputWriter::Run()
{
  GlobalData::GetInstance()->PrintNodes();
}