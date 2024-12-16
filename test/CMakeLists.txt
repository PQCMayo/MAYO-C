# SPDX-License-Identifier: Apache-2.0

if (ENABLE_PARAMS_DYNAMIC)
    add_executable(mayo_test_kat test_kat.c)
    target_link_libraries(mayo_test_kat mayo_test)

    add_executable(mayo_bench bench.c m1cycles.c)
    target_link_libraries(mayo_bench mayo)

    add_executable(mayo_test_scheme test_mayo.c)
    target_link_libraries(mayo_test_scheme mayo)
    target_include_directories(mayo_test_scheme PUBLIC ../src/common)

    add_test(MAYO_1_KAT mayo_test_kat MAYO-1)
    add_test(MAYO_2_KAT mayo_test_kat MAYO-2)
    add_test(MAYO_3_KAT mayo_test_kat MAYO-3)
    add_test(MAYO_5_KAT mayo_test_kat MAYO-5)
    add_test(MAYO_1_SELFTEST mayo_test_scheme MAYO-1)
    add_test(MAYO_2_SELFTEST mayo_test_scheme MAYO-2)
    add_test(MAYO_3_SELFTEST mayo_test_scheme MAYO-3)
    add_test(MAYO_5_SELFTEST mayo_test_scheme MAYO-5)
else()
    FOREACH(MVARIANT ${MVARIANT_S})
        string(TOLOWER ${MVARIANT} MVARIANT_LOWER)
        add_executable(mayo_test_kat_${MVARIANT} test_kat.c)
        target_link_libraries(mayo_test_kat_${MVARIANT} ${MVARIANT_LOWER}_test)

        add_executable(mayo_bench_${MVARIANT} bench.c m1cycles.c)
        target_link_libraries(mayo_bench_${MVARIANT} ${MVARIANT_LOWER})

        add_executable(mayo_test_scheme_${MVARIANT} test_mayo.c)
        target_link_libraries(mayo_test_scheme_${MVARIANT} ${MVARIANT_LOWER})
        target_include_directories(mayo_test_scheme_${MVARIANT} PUBLIC ../src/common ${INC_PLATFORM})

        add_executable(mayo_test_sample_solution_${MVARIANT} test_sample_solution.c)
        target_link_libraries(mayo_test_sample_solution_${MVARIANT} ${MVARIANT_LOWER}_test_nistapi)
        target_include_directories(mayo_test_sample_solution_${MVARIANT} PRIVATE ../include ../src/${MVARIANT_LOWER})

        add_executable(mayo_test_various_${MVARIANT} test_various.c)
        target_link_libraries(mayo_test_various_${MVARIANT} ${MVARIANT_LOWER})

        target_include_directories(mayo_test_various_${MVARIANT} PUBLIC ../src/common ../src/generic ${INC_PLATFORM})

        add_test(${MVARIANT}_KAT mayo_test_kat_${MVARIANT} ${MVARIANT})
        add_test(${MVARIANT}_SELFTEST mayo_test_scheme_${MVARIANT} ${MVARIANT})
        add_test(${MVARIANT}_VARIOUS mayo_test_various_${MVARIANT} ${MVARIANT})
        add_test(${MVARIANT}_SAMPLE_SOLUTION mayo_test_sample_solution_${MVARIANT})
    ENDFOREACH()

endif()

