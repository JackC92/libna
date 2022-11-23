cmake_minimum_required(VERSION 3.24.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_NO_SYSTEM_FROM_IMPORTED 1) # include directories using -i instead of -isystem

project("NumericalAlgorithms")

option(ENABLE_CHOLMOD "Enable CHOLMOD linear system solver" ON)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

### Compiler options
set(CMAKE_EXPORT_COMPILE_COMMANDS 1) # Emit a compile flags file to support completion engines, ignored by generators other than Makefile Generators and Ninja.

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # Using Clang or GCC
    message(STATUS "Using Clang/gcc compiler flags")
    set(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra -g3")
    set(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")
    set(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")
    
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        message(STATUS "Setting Clang specific options")
        set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
        set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fno-limit-debug-info")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        message(STATUS "Setting gcc specific options")
        set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
        set(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
    endif()
    
    set(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
    set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
    set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE        "-march=native -O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    
    add_definitions(-D_USE_MATH_DEFINES)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # Using Microsoft Visual C++
    message(STATUS "Using Microsoft Visual C++ compiler flags")
    set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # Parallel build
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
    SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
    
    add_definitions(-D_USE_MATH_DEFINES)
else()
    # Unrecognized compiler
    message(FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]")
endif()

add_subdirectory("external/googletest")

include_directories("include/")

# Add project files
file(GLOB_RECURSE HEADER_FILES FOLLOW_SYMLINKS include/*.h)
file(GLOB_RECURSE SRC_FILES FOLLOW_SYMLINKS src/*.cpp)
file(GLOB_RECURSE TEST_FILES tests/test_*.cpp)

add_library(core ${SRC_FILES} ${HEADER_FILES})
set_target_properties(core PROPERTIES LINKER_LANGUAGE CXX)

add_executable(libna "libna.cpp")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT libna)
target_link_libraries(libna PRIVATE core)

add_executable(tests "tests.cpp" ${TEST_FILES})
target_link_libraries(tests PRIVATE core gtest gtest_main)
target_include_directories(tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/include")
add_test(NAME core_tests COMMAND tests)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    foreach(header IN LISTS HEADER_FILES)
        get_filename_component(header_path "${header}" PATH)
        file(RELATIVE_PATH header_path_rel "${CMAKE_CURRENT_SOURCE_DIR}/include/na/" "${header_path}")
        string(REPLACE "/" "\\" header_path_rel_msvc "${header_path_rel}")
        source_group("Header Files\\${header_path_rel_msvc}" FILES "${header}")
    endforeach()
    foreach(source IN LISTS SRC_FILES)
        get_filename_component(source_path "${source}" PATH)
        file(RELATIVE_PATH source_path_rel "${CMAKE_CURRENT_SOURCE_DIR}/src/na/" "${source_path}")
        string(REPLACE "/" "\\" source_path_rel_msvc "${source_path_rel}")
        source_group("Source Files\\${source_path_rel_msvc}" FILES "${source}")
    endforeach()
endif()

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
target_link_libraries(core PUBLIC Eigen3::Eigen)

if (ENABLE_CHOLMOD)
    find_package(SuiteSparse CONFIG REQUIRED)
    target_link_libraries(core PUBLIC SuiteSparse::cholmod)
endif()
