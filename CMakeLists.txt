# @CMake Configuration
# Modified from aselabs
# Lab code Available at https://github.com/NCCA/aselabs-Radekeo/blob/main/Particle/CMakeLists.txt

# Try to use a version > 3.1 if avaliable
cmake_minimum_required(VERSION 3.2)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
   set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()
# Project name; it is best to use something different from the exe name
project(Fluid_build)
find_package(fmt CONFIG REQUIRED)

# Set the C++ standard to use
set(CMAKE_CXX_STANDARD 17)
# add include paths
include_directories(include)
# Add target executable and the file it is built from.
add_executable(Fluid)
target_sources(Fluid PRIVATE src/main.cpp src/Simulator.cpp src/Vec3.cpp
include/Fluid.h include/Simulator.h )
target_link_libraries(Fluid PRIVATE  fmt::fmt-header-only)

#################################################################################
# Testing code
#################################################################################

find_package(GTest CONFIG REQUIRED)
include(GoogleTest)
enable_testing()
add_executable(FluidTests)
target_sources(FluidTests PRIVATE tests/FluidTests.cpp  src/Simulator.cpp src/Vec3.cpp include/Simulator.h include/Fluid.h include/Vec3.h)
target_link_libraries(FluidTests PRIVATE GTest::gtest GTest::gtest_main )
gtest_discover_tests(FluidTests)
