# SPDX-License-Identifier: Apache-2.0

# One set of test binaries, on the libraries LIB (system randomness) and LIB_TEST
# (deterministic randomness, for the KATs).
function(mayo_test_binaries SUFFIX LIB LIB_TEST)
    add_executable(mayo_test_kat${SUFFIX} test_kat.c)
    target_link_libraries(mayo_test_kat${SUFFIX} ${LIB_TEST})

    add_executable(mayo_bench${SUFFIX} bench.c m1cycles.c)
    target_link_libraries(mayo_bench${SUFFIX} ${LIB})

    add_executable(mayo_test_scheme${SUFFIX} test_mayo.c)
    target_link_libraries(mayo_test_scheme${SUFFIX} ${LIB})
    target_include_directories(mayo_test_scheme${SUFFIX} PUBLIC ../src/common ${INC_PLATFORM})

    add_executable(mayo_test_various${SUFFIX} test_various.c)
    target_link_libraries(mayo_test_various${SUFFIX} ${LIB})
    target_include_directories(mayo_test_various${SUFFIX} PUBLIC ../src/common ../src/generic ${INC_PLATFORM})
endfunction()

# With dynamic parameters one set of binaries serves every parameter set, which is
# passed as an argument; otherwise every parameter set gets its own binaries.
if (ENABLE_PARAMS_DYNAMIC)
    mayo_test_binaries("" mayo mayo_test)
endif()

foreach(MVARIANT ${MVARIANT_S})
    string(TOLOWER ${MVARIANT} MVARIANT_LOWER)
    string(REPLACE "_" "-" MVARIANT_ARG ${MVARIANT})
    if (ENABLE_PARAMS_DYNAMIC)
        set(SUFFIX "")
    else()
        set(SUFFIX _${MVARIANT})
        mayo_test_binaries(${SUFFIX} ${MVARIANT_LOWER} ${MVARIANT_LOWER}_test)
    endif()

    add_executable(mayo_test_sample_solution_${MVARIANT} test_sample_solution.c)
    target_link_libraries(mayo_test_sample_solution_${MVARIANT} ${MVARIANT_LOWER}_test_nistapi)
    target_include_directories(mayo_test_sample_solution_${MVARIANT} PRIVATE ../include ../src/${MVARIANT_LOWER})

    add_test(${MVARIANT}_KAT mayo_test_kat${SUFFIX} ${MVARIANT_ARG})
    add_test(${MVARIANT}_SELFTEST mayo_test_scheme${SUFFIX} ${MVARIANT_ARG})
    add_test(${MVARIANT}_VARIOUS mayo_test_various${SUFFIX} ${MVARIANT_ARG})
    add_test(${MVARIANT}_SAMPLE_SOLUTION mayo_test_sample_solution_${MVARIANT})
endforeach()
