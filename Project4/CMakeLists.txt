cmake_minimum_required(VERSION 3.10)

PROJECT(proj4)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj4 proj4.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(proj4 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj4 vtkHybrid)
endif()



