find_package(Eigen3 3.4 REQUIRED CONFIG) # calling before libigl to use custom version of Eigen, otherwise libigl fetches Eigen at tags/3.4.0

add_subdirectory("external/fast_float"      EXCLUDE_FROM_ALL)
add_subdirectory("external/fmt"             EXCLUDE_FROM_ALL)
add_subdirectory("external/googletest"      EXCLUDE_FROM_ALL)
add_subdirectory("external/ImGuiFileDialog" EXCLUDE_FROM_ALL)
add_subdirectory("external/implot"          EXCLUDE_FROM_ALL)
add_subdirectory("external/libigl"          EXCLUDE_FROM_ALL)
add_subdirectory("external/matplotlib-cpp"  EXCLUDE_FROM_ALL)
add_subdirectory("external/polyscope"       EXCLUDE_FROM_ALL)

target_link_libraries(core PUBLIC Eigen3::Eigen fast_float fmt)
target_link_libraries(libna PRIVATE core igl::core imgui_filedialog implot matplotlib_cpp polyscope)
target_link_libraries(tests PRIVATE core gtest gtest_main)

set(LIBNA_HEADER_FILES "")
set(LIBNA_SOURCE_FILES "")
set(LIBNA_TEST_FILES "")
set(LIBNA_MAIN_FILE "libna.cpp")
set(LIBNA_TEST_FILE "tests.cpp")

set(_linked_targets "")

if(WITH_CHOLMOD)
    find_package(SuiteSparse_config REQUIRED CONFIG)
    find_package(CHOLMOD REQUIRED CONFIG)
    target_link_libraries(core PUBLIC SuiteSparse::CHOLMOD)
    
    list(APPEND _linked_targets "CHOLMOD")
endif()

if(WITH_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
    
        set(_CUDA_flags "")
        list(APPEND _CUDA_flags ${CMAKE_CUDA_FLAGS})
        list(APPEND _CUDA_flags "--expt-relaxed-constexpr")
        list(APPEND _CUDA_flags "--extended-lambda")
        list(JOIN _CUDA_flags " " _CUDA_flags)
        set(CMAKE_CUDA_FLAGS ${_CUDA_flags})
        
        find_package(CUDAToolkit REQUIRED)
        target_link_libraries(core PUBLIC CUDA::cublas CUDA::curand CUDA::cusolver)
        set_target_properties(core libna PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        set_target_properties(core libna PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
        
        na_include(cuda)
        set(LIBNA_MAIN_FILE "libna.cu")
        
        list(APPEND _linked_targets "CUDA")
    endif()
endif()

if(WITH_MATLAB)
    find_package(MATLAB REQUIRED)
    target_link_libraries(core PUBLIC ${Matlab_LIBRARIES})
    target_include_directories(core PUBLIC ${Matlab_INCLUDE_DIRS})
    
    list(APPEND _linked_targets "MATLAB")
endif()

if(WITH_MKL)
    # Eigen requires the LP64 interface
    set(MKL_INTERFACE "lp64")
    
    find_package(MKL REQUIRED CONFIG)
    target_link_libraries(core PUBLIC $<LINK_ONLY:MKL::MKL>)
    target_compile_options(core PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_compile_definitions(core PUBLIC "EIGEN_USE_MKL_ALL")
    target_include_directories(core PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    
    list(APPEND _linked_targets "MKL")
endif()

if(WITH_MOSEK)
    find_package(MOSEK REQUIRED)
    target_link_libraries(core PUBLIC ${MOSEK_LIBRARY})
    target_include_directories(core PUBLIC ${MOSEK_INCLUDE_DIR})
    
    na_include(mosek)
    
    list(APPEND _linked_targets "MOSEK")
endif()

if(WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(core PUBLIC $<$<CONFIG:Release>:OpenMP::OpenMP_CXX>)
    
    list(APPEND _linked_targets "OPENMP")
endif()

if(WITH_TORCH)
    find_package(Torch REQUIRED CONFIG)
    target_link_libraries(core PUBLIC torch)
    include(TorchFix)
    
    list(APPEND _linked_targets "TORCH")
endif()

if(WITH_TENSORFLOW)
    find_package(TensorFlow REQUIRED CONFIG)
    target_link_libraries(core PUBLIC TensorFlow::TensorFlow)
    
    list(APPEND _linked_targets "TENSORFLOW")
endif()

foreach(_linked_target IN LISTS _linked_targets)
    string(TOUPPER ${_linked_target} _linked_target)
    target_compile_definitions(core PUBLIC "NA_USE_${_linked_target}")
endforeach()

na_include(core)
na_include(core_test)

include(MSVCFilter)

target_sources(core PRIVATE ${LIBNA_HEADER_FILES} ${LIBNA_SOURCE_FILES})
target_sources(libna PRIVATE ${LIBNA_MAIN_FILE})
target_sources(tests PRIVATE ${LIBNA_TEST_FILE} ${LIBNA_TEST_FILES})
target_include_directories(core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
