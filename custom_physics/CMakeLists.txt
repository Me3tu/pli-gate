# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
#-------------------------------------------------------------------------------
# CMakeLists.txt for PLI Custom Physics Library
#-------------------------------------------------------------------------------
#
# Purpose:
#   Build custom Geant4 physics library for Positronium Lifetime Imaging
#   This library provides:
#     - G4PositroniumAnnihilation: Custom discrete process for Ps-based annihilation
#     - PLIPhysicsList: Complete physics list with custom annihilation
#
# Usage:
#   mkdir build
#   cd build
#   cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
#   make -j$(nproc)
#   make install
#
# Requirements:
#   - Geant4 11.x or later
#   - CMake 3.16 or later
#
#-------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

#-------------------------------------------------------------------------------
# Project Definition
#-------------------------------------------------------------------------------

project(PLICustomPhysics
    VERSION 1.0.0
    DESCRIPTION "Custom Geant4 physics for Positronium Lifetime Imaging"
    LANGUAGES CXX
)

# Require C++17 (matching Geant4 requirements)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#-------------------------------------------------------------------------------
# Find Required Packages
#-------------------------------------------------------------------------------

# Find Geant4
find_package(Geant4 REQUIRED)

if(NOT Geant4_FOUND)
    message(FATAL_ERROR "Geant4 not found! Please set Geant4_DIR or source geant4.sh")
endif()

message(STATUS "Found Geant4 ${Geant4_VERSION} at ${Geant4_DIR}")

# Include Geant4 configuration
include(${Geant4_USE_FILE})

# nlohmann/json (header-only). Required by PLIMaterialLoader.cc to parse
# materials/pli_materials.json. Header lives at /usr/include/nlohmann/
# on Arch (package nlohmann-json). find_package gives us a version check
# and the imported target nlohmann_json::nlohmann_json. If find_package
# fails (e.g. on a system without the CMake config files installed) we
# fall back to a plain include-path check.
# The bundled copy is preferred so the package builds with no system
# package installed. Set PLI_USE_SYSTEM_JSON=ON to use a system copy instead.
option(PLI_USE_SYSTEM_JSON "Use a system nlohmann/json instead of the bundled header" OFF)
if(NOT PLI_USE_SYSTEM_JSON AND EXISTS ${PROJECT_SOURCE_DIR}/third_party/nlohmann/json.hpp)
    set(NLOHMANN_JSON_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party)
    set(nlohmann_json_FOUND FALSE)
    message(STATUS "nlohmann/json: using the bundled header in third_party/")
else()
find_package(nlohmann_json 3.2.0 QUIET)
if(NOT nlohmann_json_FOUND)
    find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp
              PATHS /usr/include /usr/local/include)
    if(NOT NLOHMANN_JSON_INCLUDE_DIR)
        message(FATAL_ERROR
            "nlohmann/json not found. Install via 'pacman -S nlohmann-json' "
            "or set CMAKE_PREFIX_PATH to a directory containing it.")
    endif()
    message(STATUS "nlohmann/json found at ${NLOHMANN_JSON_INCLUDE_DIR}")
else()
    message(STATUS "nlohmann/json ${nlohmann_json_VERSION} (CMake target)")
endif()
endif()

#-------------------------------------------------------------------------------
# Include Directories
#-------------------------------------------------------------------------------

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIRS}
)

#-------------------------------------------------------------------------------
# Source Files
#-------------------------------------------------------------------------------

# List files explicitly (must match actual sources used in GATE)
set(LIBRARY_SOURCES
    ${PROJECT_SOURCE_DIR}/src/G4PositroniumAnnihilation.cc
    ${PROJECT_SOURCE_DIR}/src/G4PsAtRestProcess.cc
    ${PROJECT_SOURCE_DIR}/src/PLIPhysicsList.cc
    ${PROJECT_SOURCE_DIR}/src/PLISteppingAction.cc
    ${PROJECT_SOURCE_DIR}/src/PLIMaterialLoader.cc
    ${PROJECT_SOURCE_DIR}/src/PLIOrePowell.cc
)

set(LIBRARY_HEADERS
    ${PROJECT_SOURCE_DIR}/include/G4PositroniumAnnihilation.hh
    ${PROJECT_SOURCE_DIR}/include/G4PsAtRestProcess.hh
    ${PROJECT_SOURCE_DIR}/include/PLIPhysicsList.hh
    ${PROJECT_SOURCE_DIR}/include/PLISteppingAction.hh
    ${PROJECT_SOURCE_DIR}/include/PLIMaterialLoader.hh
    ${PROJECT_SOURCE_DIR}/include/PLIOrePowell.hh
)

# Verify all files exist
foreach(src_file ${LIBRARY_SOURCES})
    if(NOT EXISTS ${src_file})
        message(FATAL_ERROR "Source file not found: ${src_file}")
    endif()
endforeach()

foreach(hdr_file ${LIBRARY_HEADERS})
    if(NOT EXISTS ${hdr_file})
        message(FATAL_ERROR "Header file not found: ${hdr_file}")
    endif()
endforeach()

message(STATUS "Building library with ${CMAKE_BUILD_TYPE} configuration")
message(STATUS "Source files:")
foreach(src ${LIBRARY_SOURCES})
    message(STATUS "  ${src}")
endforeach()

#-------------------------------------------------------------------------------
# Library Target
#-------------------------------------------------------------------------------

# Create shared library
add_library(PLICustomPhysics SHARED
    ${LIBRARY_SOURCES}
    ${LIBRARY_HEADERS}
)

# Set library properties
set_target_properties(PLICustomPhysics PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER "${LIBRARY_HEADERS}"
    POSITION_INDEPENDENT_CODE ON
)

# Link against Geant4 libraries
target_link_libraries(PLICustomPhysics
    PUBLIC
        ${Geant4_LIBRARIES}
)

# nlohmann/json is header-only. Wire it in via either the imported
# CMake target (preferred) or a raw include path (fallback).
if(nlohmann_json_FOUND)
    target_link_libraries(PLICustomPhysics PRIVATE nlohmann_json::nlohmann_json)
else()
    target_include_directories(PLICustomPhysics PRIVATE ${NLOHMANN_JSON_INCLUDE_DIR})
endif()

# Set include directories for the target
target_include_directories(PLICustomPhysics
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

#-------------------------------------------------------------------------------
# test_load executable
#-------------------------------------------------------------------------------
#
# Standalone main() that resolves the JSON path, parses it via
# PLIMaterialLoader::loadFromJson, and prints the per-material tau_eff
# and I3 values. Run it before any GATE run: it is an end-to-end check of
# the parser, the schema validation and the Geant4 build linkage. It is NOT
# a physics check - it re-reads the same file it reports.
#
# Run with:   ./test_load
# Or:         PLI_MATERIALS_JSON=/abs/path ./test_load
#-------------------------------------------------------------------------------

add_executable(test_load tests/test_load.cc)
target_link_libraries(test_load
    PRIVATE
        PLICustomPhysics
        ${Geant4_LIBRARIES}
)
target_include_directories(test_load
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

#-------------------------------------------------------------------------------
# Compiler Warnings
#-------------------------------------------------------------------------------

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(PLICustomPhysics PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wno-unused-parameter
    )
endif()

#-------------------------------------------------------------------------------
# Installation
#-------------------------------------------------------------------------------

# Install library
install(TARGETS PLICustomPhysics
    EXPORT PLICustomPhysicsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include
)

# Install header files explicitly
install(FILES ${LIBRARY_HEADERS}
    DESTINATION include
)

# Create and install CMake package configuration files
include(CMakePackageConfigHelpers)

# Generate version file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/PLICustomPhysicsConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Generate config file
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/PLICustomPhysicsConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/PLICustomPhysicsConfig.cmake"
    INSTALL_DESTINATION lib/cmake/PLICustomPhysics
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# Install package configuration files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/PLICustomPhysicsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/PLICustomPhysicsConfigVersion.cmake"
    DESTINATION lib/cmake/PLICustomPhysics
)

# Export targets
install(EXPORT PLICustomPhysicsTargets
    FILE PLICustomPhysicsTargets.cmake
    NAMESPACE PLI::
    DESTINATION lib/cmake/PLICustomPhysics
)

#-------------------------------------------------------------------------------
# Build Summary
#-------------------------------------------------------------------------------

message(STATUS "")
message(STATUS "PLI Custom Physics Configuration Summary:")
message(STATUS "  Project Version:    ${PROJECT_VERSION}")
message(STATUS "  Geant4 Version:     ${Geant4_VERSION}")
message(STATUS "  C++ Standard:       C++${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type:         ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install Prefix:     ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  Library Name:       libPLICustomPhysics")
message(STATUS "")

#-------------------------------------------------------------------------------
# End of CMakeLists.txt
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# test_orepowell - unit test for the three-photon o-Ps final state
#-------------------------------------------------------------------------------
add_executable(test_orepowell tests/test_orepowell.cc)
target_link_libraries(test_orepowell
    PRIVATE
        PLICustomPhysics
        ${Geant4_LIBRARIES}
)
target_include_directories(test_orepowell
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)
