cmake_minimum_required(VERSION 3.10)

PROJECT(proj7)

find_package(VTK REQUIRED)
set(BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "-fopenmp")
#include(${VTK_USE_FILE})

#add_executable(astro astro.cxx)
#target_link_libraries(astro ${VTK_LIBRARIES})
add_executable(final_project final_project.cxx)
target_link_libraries(final_project ${VTK_LIBRARIES})
