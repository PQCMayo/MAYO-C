## NIST KAT generation apps
foreach(MVARIANT ${MVARIANT_S})
    string(TOLOWER ${MVARIANT} MVARIANT_LOWER)
    add_executable(PQCgenKAT_sign_${MVARIANT_LOWER} PQCgenKAT_sign.c)
    target_link_libraries(PQCgenKAT_sign_${MVARIANT_LOWER} PRIVATE ${MVARIANT_LOWER}_test_nistapi)
    target_include_directories(PQCgenKAT_sign_${MVARIANT_LOWER} PRIVATE ../include)
endforeach()

## Examples with MAYO API
if (ENABLE_PARAMS_DYNAMIC)
    add_executable(example_mayo example.c)
    target_link_libraries(example_mayo PRIVATE mayo)
    target_include_directories(example_mayo PRIVATE ../include)
else()
    foreach(MVARIANT ${MVARIANT_S})
        string(TOLOWER ${MVARIANT} MVARIANT_LOWER)
        add_executable(example_${MVARIANT_LOWER} example.c)
        target_link_libraries(example_${MVARIANT_LOWER} PRIVATE ${MVARIANT_LOWER})
        target_include_directories(example_${MVARIANT_LOWER} PRIVATE ../include)
        target_compile_definitions(example_${MVARIANT_LOWER} PRIVATE MAYO_VARIANT=${MVARIANT})
    endforeach()
endif()

## Examples with NIST API
foreach(MVARIANT ${MVARIANT_S})
    string(TOLOWER ${MVARIANT} MVARIANT_LOWER)
    add_executable(example_nistapi_${MVARIANT_LOWER} example_nistapi.c)
    target_link_libraries(example_nistapi_${MVARIANT_LOWER} PRIVATE ${MVARIANT_LOWER}_nistapi)
    target_include_directories(example_nistapi_${MVARIANT_LOWER} PRIVATE ../include ../src/${MVARIANT_LOWER})
endforeach()

