add_library(fast_float INTERFACE)
add_library(FastFloat::fast_float ALIAS fast_float)
target_include_directories(
    fast_float
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/fast_float/include)
target_compile_features(fast_float INTERFACE cxx_std_11)
