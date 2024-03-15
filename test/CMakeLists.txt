# SPDX-License-Identifier: Apache-2.0

if (ENABLE_PARAMS_DYNAMIC)
    add_executable(mayo_test_kat test_kat.c)
    target_link_libraries(mayo_test_kat mayo_test)

    add_executable(mayo_bench bench.c)
    target_link_libraries(mayo_bench mayo)

    add_executable(mayo_test_scheme test_mayo.c)
    target_link_libraries(mayo_test_scheme mayo)
    target_include_directories(mayo_test_scheme PUBLIC ../src/common)

    add_test(MAYO_1_KAT mayo_test_kat MAYO_1)
    add_test(MAYO_2_KAT mayo_test_kat MAYO_2)
    add_test(MAYO_3_KAT mayo_test_kat MAYO_3)
    add_test(MAYO_5_KAT mayo_test_kat MAYO_5)
    add_test(MAYO_1_SELFTEST mayo_test_scheme MAYO_1)
    add_test(MAYO_2_SELFTEST mayo_test_scheme MAYO_2)
    add_test(MAYO_3_SELFTEST mayo_test_scheme MAYO_3)
    add_test(MAYO_5_SELFTEST mayo_test_scheme MAYO_5)
else()
    FOREACH(MVARIANT ${MVARIANT_S})
        string(TOLOWER ${MVARIANT} MVARIANT_LOWER)
        add_executable(mayo_test_kat_${MVARIANT} test_kat.c)
        target_link_libraries(mayo_test_kat_${MVARIANT} ${MVARIANT_LOWER}_test)

        add_executable(mayo_bench_${MVARIANT} bench.c)
        target_link_libraries(mayo_bench_${MVARIANT} ${MVARIANT_LOWER})

        add_executable(mayo_test_scheme_${MVARIANT} test_mayo.c)
        target_link_libraries(mayo_test_scheme_${MVARIANT} ${MVARIANT_LOWER})
        target_include_directories(mayo_test_scheme_${MVARIANT} PUBLIC ../src/common ${INC_PLATFORM})

        add_executable(mayo_bench_table3_${MVARIANT} bench_mayo_table3.c)
        target_link_libraries(mayo_bench_table3_${MVARIANT} ${MVARIANT_LOWER})

        add_test(${MVARIANT}_KAT mayo_test_kat_${MVARIANT} ${MVARIANT})
        add_test(${MVARIANT}_SELFTEST mayo_test_scheme_${MVARIANT} ${MVARIANT})
    ENDFOREACH()

    add_executable(mayo_test_sample_solution_MAYO_1 test_sample_solution.c)
    target_link_libraries(mayo_test_sample_solution_MAYO_1 mayo_1_test_nistapi)
    target_include_directories(mayo_test_sample_solution_MAYO_1 PRIVATE ../include ../src/mayo_1)
    add_test(MAYO_1_SAMPLE_SOLUTION mayo_test_sample_solution_MAYO_1)
endif()