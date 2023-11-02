add_library(fast_float INTERFACE)
add_library(FastFloat::fast_float ALIAS fast_float)
target_include_directories(
    fast_float
    INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(fast_float INTERFACE cxx_std_11)
