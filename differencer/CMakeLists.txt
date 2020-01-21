cmake_minimum_required(VERSION 3.10)

PROJECT(differencer)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(differencer differencer.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(differencer ${VTK_LIBRARIES})
else()
  target_link_libraries(differencer vtkHybrid)
endif()



