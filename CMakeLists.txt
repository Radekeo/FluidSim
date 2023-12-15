cmake_minimum_required(VERSION 3.12)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
    set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()

# Name of the project
project(FluidSimBuild)
# This is the name of the Exe change this and it will change everywhere
set(TargetName FluidSim)
# This will include the file NGLConfig.cmake, you need to add the location to this either using
# -DCMAKE_PREFIX_PATH=~/NGL or as a system environment variable. 
find_package(NGL CONFIG REQUIRED)
# Instruct CMake to run moc automatically when needed (Qt projects only)
set(CMAKE_AUTOMOC ON)
# find Qt libs
# find Qt libs first we check for Version 6
find_package(Qt6 COMPONENTS OpenGL Widgets QUIET )
if ( Qt6_FOUND )
    message("Found Qt6 Using that")
else()
    message("Found Qt5 Using that")
    find_package(Qt5 COMPONENTS OpenGL Widgets REQUIRED)
endif()

# use C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
# Set the name of the executable we want to build
add_executable(${TargetName})

# Add NGL include path
include_directories(include $ENV{HOME}/NGL/include)
target_sources(${TargetName} PRIVATE ${PROJECT_SOURCE_DIR}/src/main.cpp
        ${PROJECT_SOURCE_DIR}/src/NGLScene.cpp
        ${PROJECT_SOURCE_DIR}/src/NGLSceneMouseControls.cpp
        ${PROJECT_SOURCE_DIR}/src/Simulator.cpp
        ${PROJECT_SOURCE_DIR}/include/Fluid.h
        ${PROJECT_SOURCE_DIR}/include/Simulator.h
        ${PROJECT_SOURCE_DIR}/include/NGLScene.h
        ${PROJECT_SOURCE_DIR}/include/WindowParams.h
)

target_link_libraries(${TargetName} PRIVATE Qt::Widgets Qt::OpenGL NGL)

add_custom_target(${TargetName}CopyShaders ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/shaders
        $<TARGET_FILE_DIR:${TargetName}>/shaders
)

##################################################################################
## Testing code
##################################################################################
#
#find_package(GTest CONFIG REQUIRED)
#include(GoogleTest)
#enable_testing()
#add_executable(FluidSimTest)
#target_sources(FluidSimTest PRIVATE tests/FluidTests.cpp
#                src/Simulator.cpp src/NGLScene.cpp src/NGLSceneMouseControls.cpp src/main.cpp
#                include/Simulator.h include/Particle.h include/NGLScene.h include/WindowParams.h
#)
#target_link_libraries(FluidSimTest PRIVATE GTest::gtest GTest::gtest_main )
#gtest_discover_tests(FluidSimTest)
