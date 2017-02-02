cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(godsim)

# GEANT4 core
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# CADMesh
find_package(cadmesh)

set(CMAKE_BUILD_TYPE release)

# User code
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(godsim godsim.cc ${sources} ${headers})
target_link_libraries(godsim ${Geant4_LIBRARIES})
target_link_libraries(godsim ${cadmesh_LIBRARIES})

install(TARGETS godsim DESTINATION ${PROJECT_SOURCE_DIR})
