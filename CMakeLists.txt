cmake_minimum_required(VERSION 3.5)

project(ParallelNonlinearFEM LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-std=c++17 -g")

set(BUILD_LIB_TYPE OBJECT)

# set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Degub)
message(STATUS "======== CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE} ========")

# --------------------------- MPI Pacakage ------------------------------
if(NOT mpi)
  set(mpi "/home/wangjianhua/Downloads/PETSc/petsc-main/test")
  message(STATUS "======== mpi PATH = ${mpi} ========")
  # set(MPI_INC_DIR "/home/wangjianhua/Downloads/PETSc/petsc-main/test/include")
  # set(MPI_LINK_DIR "/home/wangjianhua/Downloads/PETSc/petsc-main/test/lib")
endif()
include_directories(${mpi}/include/)
link_directories(${mpi}/lib/)

# --------------------------- PETSc Pacakage ------------------------------
if(NOT petsc)
  # set(petsc "${ENV}/Downloads/PETSc/petsc-main")
  set(petsc "/home/wangjianhua/Downloads/PETSc/petsc-main")
  message(STATUS "======== petsc PATH = ${petsc} ========")
endif()
# set(PTESC_INC_DIR "/home/wangjianhua/Downloads/PETSc/petsc-main/include/")
# set(PETSC_LINK_DIR "/home/wangjianhua/Downloads/PETSc/petsc-main/lib/petsc")
include_directories(${petsc}/include/)
link_directories(${petsc}/lib/petsc/)

message(STATUS "======== CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR} ========")
add_subdirectory(${CMAKE_SOURCE_DIR}/software/elements)
add_subdirectory(${CMAKE_SOURCE_DIR}/software/fem)
add_subdirectory(${CMAKE_SOURCE_DIR}/software/io)
add_subdirectory(${CMAKE_SOURCE_DIR}/software/materials)
add_subdirectory(${CMAKE_SOURCE_DIR}/software/solvers)
add_subdirectory(${CMAKE_SOURCE_DIR}/software/util)
include_directories(${CMAKE_SOURCE_DIR}/software/nlohmann)

add_library(libfemmodule SHARED
            $<TARGET_OBJECTS:fem>
            $<TARGET_OBJECTS:io>
            $<TARGET_OBJECTS:solvers>
            $<TARGET_OBJECTS:util>)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME}
                      mpi
                      petsc
                      libfemmodule
                      materials
                      elements)

message(STATUS "======== BUILDING " ${PROJECT_NAME} " SUCCESSFUL! ========")
