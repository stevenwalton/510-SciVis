cmake_minimum_required(VERSION 3.10)

PROJECT(proj6)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(BUILD_TYPE "Debug")

add_executable(proj6 proj6.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(proj6 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj6 vtkHybrid)
endif()



