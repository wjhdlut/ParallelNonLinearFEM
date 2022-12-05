
#include <io/MeshWriter.h>
#include <util/DataStructure.h>

MeshWriter::MeshWriter(const nlohmann::json &props) : BaseModule(props)
{
  
}

MeshWriter::~MeshWriter()
{}

void MeshWriter::Run()
{
  if(0 != (GlobalData::GetInstance()->m_cycle%m_interval)) return;

}