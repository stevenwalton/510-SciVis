cmake_minimum_required(VERSION 3.10)

PROJECT(proj2)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj2 proj2.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(proj2 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj2 vtkHybrid)
endif()



