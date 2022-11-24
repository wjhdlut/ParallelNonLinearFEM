#ifndef OBJECTFACTORY_H
#define OBJECTFACTORY_H

#include <string>
#include <unordered_map>
#include <memory>

template<class MyClass, typename ...Argtype>
void* __createObjFunc(Argtype... arg)
{
  return new MyClass(arg...);
}

#ifndef REFLECT_REGISTER
#define ReflectRegister(MyClass, ...)\
    static int __type##MyClass = ObjectFactory::registerFunction(#MyClass, \
    (void*)&__createObjFunc<MyClass, ##__VA_ARGS__>);
#endif

class ObjectFactory
{
public:
  template<class _tp, typename ..._args>
  static std::shared_ptr<_tp> CreateObject(const std::string &name, _args... args)
  {
    typedef _tp*(*funcPtr)(_args...);
    auto &func_map = getClassFunctionMap();
    auto find_func = func_map.find(name);
    if(find_func == func_map.end()) return nullptr;
    else
    {
      _tp *class_ptr = reinterpret_cast<funcPtr>(find_func->second)(args...);
      return std::shared_ptr<_tp>(class_ptr);
    }
  }

  static int registerFunction(const std::string &class_name, void *function_ptr)
  {
    getClassFunctionMap()[class_name] = function_ptr;
    return 0;
  }

private:
  static std::unordered_map<std::string, void*> &getClassFunctionMap()
  {
    static std::unordered_map<std::string, void*> class_function_map;
    return class_function_map;
  }
};

#endif // OBJECTFACTORY_H