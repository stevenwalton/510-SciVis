cmake_minimum_required(VERSION 3.10)

PROJECT(proj3)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj3 proj3.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(proj3 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj3 vtkHybrid)
endif()



