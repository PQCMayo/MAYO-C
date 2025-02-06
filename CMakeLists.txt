# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.10)
project(MAYO VERSION 1.0 LANGUAGES C CXX ASM)

set(MAYO_SO_VERSION "0")
set(CMAKE_C_STANDARD 99)

set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1")
include(CTest)

option(ENABLE_STRICT "Build with strict compile options." ON)

if(ENABLE_STRICT)
   message("Enable strict flag ON")
endif()

option(ENABLE_TESTS  "Enable compilation of tests." ON)
option(ENABLE_CT_TESTING  "Enable compilation for constant time testing." OFF)
option(SOURCE_PATH "Where the source code is stored, default in the building tree" OFF)
#option(ENABLE_PARAMS_DYNAMIC "Enable dynamic use of MAYO parameters" OFF)

SET(MVARIANT_S "MAYO_1;MAYO_2;MAYO_3;MAYO_5")

include(.cmake/flags.cmake)
include(.cmake/sanitizers.cmake)
include(.cmake/target.cmake)

add_subdirectory(src)
add_subdirectory(apps)

if(ENABLE_TESTS)
	enable_testing()
	add_subdirectory(test)
endif()
