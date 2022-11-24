
#include "nlohmann/json.hpp"

class Properties
{
public:
  Properties();
  ~Properties();

private:
  nlohmann::json m_pro;
};