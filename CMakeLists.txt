cmake_minimum_required(VERSION 3.6)

project(c-urchin
    VERSION 1.0
    DESCRIPTION "C++ Cyclotron Radiation Power Calculation in a Waveguide"
    LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE PATH "Install top-level directory" FORCE)
endif()

find_package(Boost 1.56 REQUIRED COMPONENTS filesystem program_options)

set(DATA_DIR "${PROJECT_SOURCE_DIR}/Data")
message("Data Directory: ${DATA_DIR}")
configure_file("Source/Param.hh.in" "Source/Param.hh")

# Includes the new header
include_directories(${CMAKE_BINARY_DIR}/Source)

add_executable(c-urchin
    Source/Beta.hh
    Source/Beta.cc
    Source/CircularWaveguide.hh
    Source/CircularWaveguide.cc
    Source/Constants.hh
    Source/Param.hh.in
    Source/RectangularWaveguide.hh
    Source/RectangularWaveguide.cc
    Source/Waveguide.hh
    Source/Waveguide.cc
    Source/main.cc)

target_link_libraries(c-urchin Boost::filesystem Boost::program_options)

install(TARGETS c-urchin)
