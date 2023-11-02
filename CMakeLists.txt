cmake_minimum_required(VERSION 3.24.0)

set(CMAKE_CXX_STANDARD  17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES "native")
set(CMAKE_NO_SYSTEM_FROM_IMPORTED 1) # include directories using -i instead of -isystem

project("NumericalAlgorithms" LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_LIST_DIR}/cmake"
    "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

option(WITH_CHOLMOD    "Build with CHOLMOD"    OFF)
option(WITH_CUDA       "Build with CUDA"       OFF)
option(WITH_MATLAB     "Build with MATLAB"     OFF)
option(WITH_MKL        "Build with Intel MKL"   ON)
option(WITH_MOSEK      "Build with Mosek"      OFF)
option(WITH_OPENMP     "Build with OpenMP"      ON)
option(WITH_TORCH      "Build with PyTorch"    OFF)
option(WITH_TENSORFLOW "Build with TensorFlow" OFF)

include(CheckLanguage)
check_language(CUDA)

set(_CUDA_flags "")
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

file(GLOB_RECURSE HEADER_FILES FOLLOW_SYMLINKS include/na/*.h)
file(GLOB_RECURSE SOURCE_FILES FOLLOW_SYMLINKS src/na/*.c src/na/*.cpp src/na/*.f)
file(GLOB_RECURSE TEST_FILES tests.cpp tests/na/*.h tests/na/*_test.cpp)
file(GLOB_RECURSE DRIVER_FILES include/driver/*.h src/driver/*.c src/driver/*.cpp src/driver/*.f)
file(GLOB MAIN_FILE libna.cpp)

if (CMAKE_CUDA_COMPILER AND WITH_CUDA)
    file(GLOB_RECURSE CUDA_HEADER_FILES FOLLOW_SYMLINKS include/na/*.cuh)
    file(GLOB_RECURSE CUDA_SOURCE_FILES FOLLOW_SYMLINKS src/na/*.cu)
    file(GLOB_RECURSE CUDA_TEST_FILES tests/na/*.cuh tests/na/*_test.cu)
    file(GLOB_RECURSE CUDA_DRIVER_FILES include/driver/*.cuh src/driver/*.cu)
	file(GLOB MAIN_FILE libna.cu)
	
	list(APPEND HEADER_FILES ${CUDA_HEADER_FILES})
	list(APPEND SOURCE_FILES ${CUDA_SOURCE_FILES})
	list(APPEND TEST_FILES ${CUDA_TEST_FILES})
	list(APPEND DRIVER_FILES ${CUDA_DRIVER_FILES})
endif()

include(MSVCFilter)
include(CompilerOptions)

add_library(core ${HEADER_FILES} ${SOURCE_FILES})
add_executable(libna ${MAIN_FILE} ${DRIVER_FILES})
add_executable(tests ${TEST_FILES})

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT libna)
set_target_properties(core libna PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(core libna PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(core libna PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googlemock/include")
target_include_directories(tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/include")

find_package(Eigen3 3.4 REQUIRED CONFIG) # calling before libigl to use custom version of Eigen, otherwise libigl fetches Eigen at tags/3.4.0
find_package(TBB REQUIRED CONFIG)

add_subdirectory("external/fast_float"      EXCLUDE_FROM_ALL)
add_subdirectory("external/fmt"             EXCLUDE_FROM_ALL)
add_subdirectory("external/googletest"      EXCLUDE_FROM_ALL)
add_subdirectory("external/ImGuiFileDialog" EXCLUDE_FROM_ALL)
add_subdirectory("external/implot"          EXCLUDE_FROM_ALL)
add_subdirectory("external/libigl"          EXCLUDE_FROM_ALL)
add_subdirectory("external/matplotlib-cpp"  EXCLUDE_FROM_ALL)
add_subdirectory("external/polyscope"       EXCLUDE_FROM_ALL)

target_link_libraries(core PUBLIC Eigen3::Eigen fast_float fmt TBB::tbb)
target_link_libraries(libna PRIVATE core igl::core imgui_filedialog implot matplotlib_cpp polyscope)
target_link_libraries(tests PRIVATE core gtest gtest_main)

set(_linked_targets "")

if (WITH_CHOLMOD)
    find_package(SuiteSparse REQUIRED CONFIG)
    target_link_libraries(core PUBLIC SuiteSparse::cholmod)
    
    list(APPEND _linked_targets "CHOLMOD")
endif()

if (WITH_CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(core PUBLIC CUDA::cublas CUDA::curand CUDA::cusolver)
    
    list(APPEND _linked_targets "CUDA")
endif()

if (WITH_MATLAB)
    find_package(MATLAB REQUIRED)
    target_link_libraries(core PUBLIC ${Matlab_LIBRARIES})
	target_include_directories(core PUBLIC ${Matlab_INCLUDE_DIRS})
    
    list(APPEND _linked_targets "MATLAB")
endif()

if (WITH_MKL)
    # Eigen requires the LP64 interface
    set(MKL_INTERFACE "lp64")
    
    find_package(MKL REQUIRED CONFIG)
    target_link_libraries(core PUBLIC $<LINK_ONLY:MKL::MKL>)
    target_compile_options(core PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
	target_compile_definitions(core PUBLIC "EIGEN_USE_MKL_ALL")
    target_include_directories(core PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    
    list(APPEND _linked_targets "MKL")
endif()

if (WITH_MOSEK)
    find_package(MOSEK REQUIRED)
	target_link_libraries(core PUBLIC ${MOSEK_LIBRARY})
	target_include_directories(core PUBLIC ${MOSEK_INCLUDE_DIR})
	
	list(APPEND _linked_targets "MOSEK")
endif()

if (WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(core PUBLIC $<$<CONFIG:Release>:OpenMP::OpenMP_CXX>)
	
	list(APPEND _linked_targets "OPENMP")
endif()

if (WITH_TORCH)
    find_package(Torch REQUIRED CONFIG)
    target_link_libraries(core PUBLIC torch)
    include(TorchFix)
    
    list(APPEND _linked_targets "TORCH")
endif()

if (WITH_TENSORFLOW)
    find_package(TensorFlow REQUIRED CONFIG)
    target_link_libraries(core PUBLIC TensorFlow::TensorFlow)
    
    list(APPEND _linked_targets "TENSORFLOW")
endif()

foreach(_linked_target IN LISTS _linked_targets)
    string(TOUPPER ${_linked_target} _linked_target)
	target_compile_definitions(core PUBLIC "NA_USE_${_linked_target}")
endforeach()
