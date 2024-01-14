SET(MODULE_NAME elements)

add_subdirectory(shapefunctions)

# find all head and source files in module elements
file(GLOB_RECURSE ${MODULE_NAME}_SRC_FILES source/*.cpp)
file(GLOB_RECURSE ${MODULE_NAME}_INC_FILES include/*.h)

# building a dynamic library
add_library(${MODULE_NAME} SHARED
            ${${MODULE_NAME}_SRC_FILES}
            $<TARGET_OBJECTS:shapefunctions>)

# reduce compile dependece of dynamic libraries
if(NOT BUILD_ALL_SHARED)
  set_target_properties(${MODULE_NAME}
                        PROPERTIES
                        POSITION_INDEPENDENT_CODE 1)
endif()