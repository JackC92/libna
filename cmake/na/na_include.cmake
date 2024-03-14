function(na_include module)
    string(TOUPPER ${module} module_uc)
    
    if(LIBNA_${module_uc}_INCLUDED)
        return()
    endif()
    
    set(MODULE_HEADER_FILES "")
    set(MODULE_SOURCE_FILES "")
    include(${NumericalAlgorithms_SOURCE_DIR}/cmake/na/modules/${module}.cmake)
    if("${module_uc}" MATCHES "TEST$")
        list(APPEND LIBNA_TEST_FILES ${MODULE_HEADER_FILES})
        list(APPEND LIBNA_TEST_FILES ${MODULE_SOURCE_FILES})
        set(LIBNA_TEST_FILES ${LIBNA_TEST_FILES} PARENT_SCOPE)
    else()
        list(APPEND LIBNA_HEADER_FILES ${MODULE_HEADER_FILES})
        list(APPEND LIBNA_SOURCE_FILES ${MODULE_SOURCE_FILES})
        set(LIBNA_HEADER_FILES ${LIBNA_HEADER_FILES} PARENT_SCOPE)
        set(LIBNA_SOURCE_FILES ${LIBNA_SOURCE_FILES} PARENT_SCOPE)
    endif()
    set(LIBNA_${module_uc}_INCLUDED ON PARENT_SCOPE)
endfunction()