set(MODULE_DIRS
    .
    core
    geometry
    graphics
    integrate
    io
    linalg
    optimize
    special
    spline
    type_traits)

foreach (dir IN LISTS MODULE_DIRS)
    file(GLOB DIR_HEADER_FILES "${NumericalAlgorithms_SOURCE_DIR}/include/na/${dir}/*.h")
    file(GLOB DIR_SOURCE_FILES 
        "${NumericalAlgorithms_SOURCE_DIR}/src/na/${dir}/*.c"
        "${NumericalAlgorithms_SOURCE_DIR}/src/na/${dir}/*.cpp")
    
    list(APPEND MODULE_HEADER_FILES ${DIR_HEADER_FILES})
    list(APPEND MODULE_SOURCE_FILES ${DIR_SOURCE_FILES})
endforeach()
