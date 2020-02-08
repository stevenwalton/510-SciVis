cmake_minimum_required(VERSION 3.10)

PROJECT(proj5)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(BUILD_TYPE "Debug")

add_executable(proj5 proj5.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(proj5 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj5 vtkHybrid)
endif()



