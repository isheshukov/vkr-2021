project(Tropical)
cmake_minimum_required(VERSION 3.15)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/lib/eigen-3.3.8/cmake)
set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/eigen-3.3.8)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(Eigen3 REQUIRED)


add_subdirectory(lib)
add_subdirectory(app)
