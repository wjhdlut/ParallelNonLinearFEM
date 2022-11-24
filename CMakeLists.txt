cmake_minimum_required(VERSION 3.5)

project(ParallelNonlinearFEM LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-std=c++17 -g")

# --------------------------- MPI Pacakage ------------------------------
set(MPI_INC_DIR "/home/wangjianhua/Downloads/PETSc/petsc-main/test/include")
set(MPI_LINK_DIR "/home/wangjianhua/Downloads/PETSc/petsc-main/test/lib")
include_directories(${MPI_INC_DIR})
link_directories(${MPI_LINK_DIR})


# --------------------------- PETSc Pacakage ------------------------------
set(PTESC_INC_DIR "/home/wangjianhua/Downloads/PETSc/petsc-main/include/")
set(PETSC_LINK_DIR "/home/wangjianhua/Downloads/PETSc/petsc-main/lib/petsc")
include_directories(${PTESC_INC_DIR})
link_directories(${PETSC_LINK_DIR})

# ----------------------------- head file ------------------------------
include_directories(${CMAKE_SOURCE_DIR}/include)

FILE(GLOB_RECURSE DIR_SRCS ${CMAKE_SOURCE_DIR}/source/*.cpp)

add_executable(${PROJECT_NAME} main.cpp ${DIR_SRCS})

target_link_libraries(${PROJECT_NAME} libmpi.so)
target_link_libraries(${PROJECT_NAME} libpetsc.so)
