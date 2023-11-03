cmake_minimum_required(VERSION 3.24.0)

set(CMAKE_CXX_STANDARD  17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES "native")
set(CMAKE_NO_SYSTEM_FROM_IMPORTED 1) # include directories using -i instead of -isystem

project("NumericalAlgorithms" LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_LIST_DIR}/cmake"
    "${CMAKE_CURRENT_LIST_DIR}/cmake/find_modules"
    "${CMAKE_CURRENT_LIST_DIR}/cmake/misc"
    "${CMAKE_CURRENT_LIST_DIR}/cmake/na"
    "${CMAKE_CURRENT_LIST_DIR}/cmake/na/modules")

option(WITH_CHOLMOD    "Build with CHOLMOD"    OFF)
option(WITH_CUDA       "Build with CUDA"       OFF)
option(WITH_MATLAB     "Build with MATLAB"     OFF)
option(WITH_MKL        "Build with Intel MKL"   ON)
option(WITH_MOSEK      "Build with Mosek"      OFF)
option(WITH_OPENMP     "Build with OpenMP"      ON)
option(WITH_TORCH      "Build with PyTorch"    OFF)
option(WITH_TENSORFLOW "Build with TensorFlow" OFF)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(CompilerOptions)

add_library(core)
add_executable(libna)

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT libna)
set_target_properties(core libna PROPERTIES LINKER_LANGUAGE CXX)

include(na_include)

include(libna)
