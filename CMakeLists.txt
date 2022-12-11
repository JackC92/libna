cmake_minimum_required(VERSION 3.24.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_NO_SYSTEM_FROM_IMPORTED 1) # include directories using -i instead of -isystem

project("NumericalAlgorithms" LANGUAGES CXX Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

option(WITH_CHOLMOD "Build with CHOLMOD" ON)
option(WITH_MATLAB "Build with MATLAB" ON)
option(WITH_MKL "Build with Intel MKL" ON)
option(WITH_OPENMP "Build with OpenMP" ON)
option(WITH_TORCH "Build with PyTorch" ON)
option(WITH_TENSORFLOW "Build with TensorFlow" ON)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory("external/googletest")

file(GLOB_RECURSE HEADER_FILES FOLLOW_SYMLINKS include/*.h)
file(GLOB_RECURSE SRC_FILES FOLLOW_SYMLINKS src/*.cpp)
file(GLOB_RECURSE TEST_FILES tests/test_*.cpp)

include(MSVCFilter)
include(CompilerOptions)

add_library(core ${HEADER_FILES} ${SRC_FILES})
add_executable(libna "libna.cpp")
add_executable(tests "tests.cpp" ${TEST_FILES})

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT libna)
set_target_properties(core PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(libna PUBLIC core)
target_link_libraries(tests PUBLIC core gtest gtest_main)

target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/include")

find_package(Eigen3 3.4 REQUIRED CONFIG)
target_link_libraries(core PUBLIC Eigen3::Eigen)

set(_linked_targets "")
if (WITH_CHOLMOD)
    find_package(SuiteSparse REQUIRED CONFIG)
	target_link_libraries(core PUBLIC SuiteSparse::cholmod)
	
	list(APPEND _linked_targets "CHOLMOD")
endif()

if (WITH_MATLAB)
	find_package(Matlab REQUIRED CONFIG)
	target_link_libraries(core PUBLIC Matlab::matlab)
	
	list(APPEND _linked_targets "MATLAB")
endif()

if (WITH_MKL)
    # Eigen requires the LP64 interface
    set(MKL_INTERFACE_FULL "intel_lp64")
	add_definitions(-DEIGEN_USE_MKL_ALL)
	
	find_package(MKL REQUIRED CONFIG)
	target_compile_options(core PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
	target_include_directories(core PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
	target_link_libraries(core PUBLIC $<LINK_ONLY:MKL::MKL>)
	
	list(APPEND _linked_targets "MKL")
endif()

if (WITH_OPENMP)
    find_package(OpenMP REQUIRED)
	target_link_libraries(core PUBLIC $<$<CONFIG:Release>:OpenMP::OpenMP_CXX>)
endif()

if (WITH_TORCH)
	find_package(Torch REQUIRED CONFIG)
	target_link_libraries(core PUBLIC torch)
	include(TorchFix)
	
	list(APPEND _linked_targets "TORCH")
endif()

if (WITH_TENSORFLOW)
	find_package(TensorFlow REQUIRED CONFIG)
	target_link_libraries(core PUBLIC TensorFlow::tensorflow)
	
	list(APPEND _linked_targets "TENSORFLOW")
endif()

foreach(_linked_target IN LISTS _linked_targets)
	add_definitions("-DNA_USE_${_linked_target}")
endforeach()
