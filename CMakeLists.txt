# @CMake Configuration
# Modified from aselabs
# Lab code Available at https://github.com/NCCA/aselabs-Radekeo/blob/main/Particle/CMakeLists.txt

# Try to use a version > 3.1 if avaliable
cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND ENV{CMAKE_TOOLCHAIN_FILE})
    set (CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()
# Project name; it is best to use something different from the exe name
project(Fluid_build)
set(CMAKE_CXX_STANDARD 17) # Set the C++ standard to use
# add include paths

# Add target executable and the file it is built from.
add_executable(Fluid
        include/Domain.h)
target_include_directories(Fluid PRIVATE include)
target_sources(Fluid PRIVATE src/main.cpp src/Simulator.cpp src/Vec3.cpp include/Fluid.h include/Simulator.h)

## OpenGL
#find_package(OpenGL REQUIRED)
#find_package(glfw3 REQUIRED)
##target_include_directories(Fluid PRIVATE ${GLFW_INCLUDE_DIRS})
#target_link_libraries(Fluid PRIVATE glfw)
##target_link_libraries(Fluid PRIVATE ${GLFW_LIBRARIES})
#
include_directories(${OPENGL_INCLUDE_DIR})

#################################################################################
# Testing code
#################################################################################

include(FetchContent)
FetchContent_Declare(
        googletest
        # Specify the commit you depend on and update it regularly.
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(FluidTests)
target_include_directories(FluidTests PRIVATE include)
target_sources(FluidTests PRIVATE tests/FluidTests.cpp  src/Simulator.cpp src/Vec3.cpp include/Simulator.h include/Fluid.h include/Vec3.h)
target_link_libraries(FluidTests gtest_main)
add_test(NAME FluidTests_t COMMAND FluidTests)
