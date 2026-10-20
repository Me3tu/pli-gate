# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
cmake_minimum_required(VERSION 3.16)
project(pli_inflight_test CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
# The SHIPPED PLI sources (sha256 checked against ENVIRONMENT_MANIFEST.json by run_tests.sh)
# v1.1.0: relative to this folder. The variant source is the shipped PLIPhysicsList.cc
# with ../gate-inflight-stock.patch applied, made at configure time in the build tree.
set(PLI ${CMAKE_CURRENT_SOURCE_DIR}/../../../custom_physics)
set(PLI_VARIANT_SRC ${CMAKE_CURRENT_BINARY_DIR}/PLIPhysicsList_inflight.cc)
configure_file(${PLI}/src/PLIPhysicsList.cc ${PLI_VARIANT_SRC} COPYONLY)
execute_process(COMMAND patch --silent ${PLI_VARIANT_SRC} ${CMAKE_CURRENT_SOURCE_DIR}/../gate-inflight-stock.patch
                RESULT_VARIABLE PLI_PATCH_RC)
if(NOT PLI_PATCH_RC EQUAL 0)
  message(FATAL_ERROR "gate-inflight-stock.patch does not apply to ${PLI}/src/PLIPhysicsList.cc")
endif()
add_executable(pli_inflight_test pli_inflight_test.cc
  ${PLI}/src/G4PositroniumAnnihilation.cc
  ${PLI}/src/G4PsAtRestProcess.cc
  ${PLI_VARIANT_SRC}
  ${PLI}/src/PLIMaterialLoader.cc
  ${PLI}/src/PLIOrePowell.cc)
target_include_directories(pli_inflight_test PRIVATE ${PLI}/include ${PLI}/third_party)
target_link_libraries(pli_inflight_test ${Geant4_LIBRARIES})
