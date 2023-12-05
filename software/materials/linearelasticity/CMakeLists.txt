SET(MODULE_NAME linearelasticity)

# find all head and source files in module elements
file(GLOB_RECURSE ${MODULE_NAME}_SRC_FILES source/*.cpp)
file(GLOB_RECURSE ${MODULE_NAME}_INC_FILES include/*.h)

# building a dynamic library
add_library(${MODULE_NAME} ${BUILD_LIB_TYPE} ${${MODULE_NAME}_SRC_FILES})

# reduce compile independece of dynamic libraries
if(NOT BUILD_ALL_SHARED)
  set_target_properties(${MODULE_NAME}
                        PROPERTIES
                        POSITION_INDEPENDENT_CODE 1)
endif()