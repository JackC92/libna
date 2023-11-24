set(MODULE_DIRS
    .)

foreach (dir IN LISTS MODULE_DIRS)
    file(GLOB DIR_HEADER_FILES "${NumericalAlgorithms_SOURCE_DIR}/tests/${dir}/*.h")
    file(GLOB DIR_SOURCE_FILES 
        "${NumericalAlgorithms_SOURCE_DIR}/tests/${dir}/*.c"
        "${NumericalAlgorithms_SOURCE_DIR}/tests/${dir}/*.cpp")
    
    list(APPEND MODULE_HEADER_FILES ${DIR_HEADER_FILES})
    list(APPEND MODULE_SOURCE_FILES ${DIR_SOURCE_FILES})
endforeach()
