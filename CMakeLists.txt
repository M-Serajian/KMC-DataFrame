cmake_minimum_required(VERSION 3.17)
project(kmcpy VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ZLIB REQUIRED)
find_package(pybind11 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(Threads REQUIRED)

find_program(MAKE_CMD NAMES make gmake REQUIRED)

set(LIBKMC_CORE "${CMAKE_CURRENT_SOURCE_DIR}/bin/libkmc_core.a")
set(KMC_RUNNER_H "${CMAKE_CURRENT_SOURCE_DIR}/include/kmc_runner.h")

# ---------------------------------------------------------------------------
# Build libkmc_core.a using the Makefile.
# Override STATIC_CFLAGS/STATIC_LFLAGS to remove -static (incompatible with
# shared library) and add -fPIC (required for linking into _core.so).
# We pass these as make variable overrides on the command line.
# ---------------------------------------------------------------------------
add_custom_command(
    OUTPUT "${LIBKMC_CORE}" "${KMC_RUNNER_H}"
    COMMAND ${MAKE_CMD}
        -j4
        bin/libkmc_core.a
        STATIC_CFLAGS=-pthread
        STATIC_LFLAGS=-pthread
        "CFLAGS=-Wall -O3 -fsigned-char -fPIC -m64 -pthread -std=c++14"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Building libkmc_core.a via Makefile (with -fPIC)"
)

add_custom_target(kmc_core_lib DEPENDS "${LIBKMC_CORE}")

pybind11_add_module(_core MODULE kmcpy/_core.cpp)

add_dependencies(_core kmc_core_lib)

target_include_directories(_core PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/kmc_core"
    "${CMAKE_CURRENT_SOURCE_DIR}/kmc_api"
    "${ZLIB_INCLUDE_DIRS}"
)

target_link_libraries(_core PRIVATE
    "${LIBKMC_CORE}"
    ZLIB::ZLIB
    Threads::Threads
)

install(TARGETS _core DESTINATION kmcpy)