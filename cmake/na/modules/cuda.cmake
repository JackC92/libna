set(MODULE_DIRS
    cuda)

foreach (dir IN LISTS MODULE_DIRS)
    file(GLOB DIR_HEADER_FILES "${NumericalAlgorithms_SOURCE_DIR}/include/na/${dir}/*.cuh")
    file(GLOB DIR_SOURCE_FILES 
        "${NumericalAlgorithms_SOURCE_DIR}/src/na/${dir}/*.cu")
    
    list(APPEND MODULE_HEADER_FILES ${DIR_HEADER_FILES})
    list(APPEND MODULE_SOURCE_FILES ${DIR_SOURCE_FILES})
endforeach()
