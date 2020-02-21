cmake_minimum_required(VERSION 3.10)

PROJECT(proj7)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj7 proj7.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(proj7 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj7 vtkHybrid)
endif()



