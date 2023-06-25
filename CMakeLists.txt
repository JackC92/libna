cmake_minimum_required(VERSION 3.24.0)

set(CMAKE_CXX_STANDARD  17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES "native")
set(CMAKE_NO_SYSTEM_FROM_IMPORTED 1) # include directories using -i instead of -isystem

project("NumericalAlgorithms" LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

option(WITH_CHOLMOD "Build with CHOLMOD" OFF)
option(WITH_CUDA "Build with CUDA" OFF)
option(WITH_MATLAB "Build with MATLAB" OFF)
option(WITH_MKL "Build with Intel MKL" ON)
option(WITH_OPENMP "Build with OpenMP" ON)
option(WITH_TORCH "Build with PyTorch" OFF)
option(WITH_TENSORFLOW "Build with TensorFlow" OFF)

set(_CUDA_flags "")
set(_linked_targets "")

include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER AND WITH_CUDA)
    enable_language(CUDA)
    
    list(APPEND _CUDA_flags ${CMAKE_CUDA_FLAGS})
    list(APPEND _CUDA_flags "--expt-relaxed-constexpr")
    list(APPEND _CUDA_flags "--extended-lambda")
    list(JOIN _CUDA_flags " " _CUDA_flags)
    set(CMAKE_CUDA_FLAGS ${_CUDA_flags})
endif()

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory("external/fmt"             EXCLUDE_FROM_ALL)
add_subdirectory("external/googletest"      EXCLUDE_FROM_ALL)
add_subdirectory("external/ImGuiFileDialog" EXCLUDE_FROM_ALL)
add_subdirectory("external/implot"          EXCLUDE_FROM_ALL)
add_subdirectory("external/libigl"          EXCLUDE_FROM_ALL)
add_subdirectory("external/matplotlib-cpp"  EXCLUDE_FROM_ALL)
add_subdirectory("external/polyscope"       EXCLUDE_FROM_ALL)

if (CMAKE_CUDA_COMPILER AND WITH_CUDA)
    file(GLOB_RECURSE HEADER_FILES FOLLOW_SYMLINKS include/na/*.h include/na/*.cuh)
    file(GLOB_RECURSE SRC_FILES FOLLOW_SYMLINKS src/na/*.c src/na/*.cpp src/na/*.cu src/.na/*f)
    file(GLOB_RECURSE TEST_FILES tests.cpp tests/na/*.h tests/na/*.cuh tests/na/*_test.cpp tests/na/*_test.cu)
    file(GLOB DRIVER_FILE libna.cpp libna.cu)
else()
    file(GLOB_RECURSE HEADER_FILES FOLLOW_SYMLINKS include/na/*.h)
    file(GLOB_RECURSE SRC_FILES FOLLOW_SYMLINKS src/na/*.c src/na/*.cpp src/na/*.f)
    file(GLOB_RECURSE TEST_FILES tests.cpp tests/na/*.h tests/na/*_test.cpp)
    file(GLOB DRIVER_FILE libna.cpp)
endif()

include(MSVCFilter)
include(CompilerOptions)

add_library(core ${HEADER_FILES} ${SRC_FILES})
add_executable(libna ${DRIVER_FILE})
add_executable(tests ${TEST_FILES})

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT libna)
set_target_properties(core PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(core libna PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(core libna PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

target_link_libraries(libna PRIVATE core fmt::fmt igl::core imgui_filedialog implot matplotlib_cpp polyscope)
target_link_libraries(tests PRIVATE core gtest gtest_main)

target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googlemock/include")
target_include_directories(tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/include")

find_package(Eigen3 3.4 REQUIRED CONFIG)
target_link_libraries(core PUBLIC Eigen3::Eigen)

find_package(TBB REQUIRED CONFIG)
target_link_libraries(core PUBLIC TBB::tbb)

if (WITH_CHOLMOD)
    find_package(SuiteSparse REQUIRED CONFIG)
    target_link_libraries(core PUBLIC SuiteSparse::cholmod)
    
    list(APPEND _linked_targets "CHOLMOD")
endif()

if (WITH_CUDA)
    find_package(CUDAToolkit)
    target_link_libraries(core PUBLIC CUDA::curand)
    
    list(APPEND _linked_targets "CUDA")
endif()

if (WITH_MATLAB)
    find_package(MATLAB REQUIRED CONFIG)
    target_link_libraries(core PUBLIC MATLAB::matlab)
    
    list(APPEND _linked_targets "MATLAB")
endif()

if (WITH_MKL)
    # Eigen requires the LP64 interface
    set(MKL_INTERFACE "lp64")
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
