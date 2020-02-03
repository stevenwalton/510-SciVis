cmake_minimum_required(VERSION 3.10)

PROJECT(proj4G)

SET(CMAKE_BUILD_TYPE "Debug")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj4G proj4G.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(proj4G ${VTK_LIBRARIES})
else()
  target_link_libraries(proj4G vtkHybrid)
endif()



