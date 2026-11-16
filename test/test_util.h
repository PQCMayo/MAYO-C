// SPDX-License-Identifier: Apache-2.0

#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include <mayo.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef ENABLE_CT_TESTING
#include <valgrind/memcheck.h>
#endif

/* Print len bytes in hex, in numbered rows of row_len bytes, or on one line if row_len
   is 0. Under constant-time testing a copy declared public is printed, so the bytes
   themselves stay tracked as secret. */
static inline void print_hex(const unsigned char *hex, int len, int row_len) {
#ifdef ENABLE_CT_TESTING
    unsigned char *copy = calloc(len, 1);
    memcpy(copy, hex, len);
    VALGRIND_MAKE_MEM_DEFINED(copy, len);
    hex = copy;
#endif
    if (row_len > 0) {
        printf("%4d:", 0);
    }
    for (int i = 0; i < len; ++i) {
        if (row_len > 0 && i % row_len == 0 && i > 0) {
            printf("\n%4d:", i / row_len);
        }
        printf("%02x", hex[i]);
    }
    printf("\n");
#ifdef ENABLE_CT_TESTING
    free(copy);
#endif
}

/* The parameter set a test program runs on: named by argv[1] ("MAYO-1", ...) with
   dynamic parameters, fixed at compile time (NULL) otherwise. Returns 0 on success. */
static inline int mayo_params_from_args(int argc, char *argv[], const mayo_params_t **p) {
#ifdef ENABLE_PARAMS_DYNAMIC
    static const struct {
        const char *name;
        const mayo_params_t *params;
    } sets[] = {
        { "MAYO-1", &MAYO_1 }, { "MAYO-2", &MAYO_2 },
        { "MAYO-3", &MAYO_3 }, { "MAYO-5", &MAYO_5 },
    };
    if (argc < 2) {
        printf("usage: %s MAYO-1|MAYO-2|MAYO-3|MAYO-5\n", argv[0]);
        return -1;
    }
    for (size_t i = 0; i < sizeof(sets) / sizeof(sets[0]); i++) {
        if (!strcmp(argv[1], sets[i].name)) {
            *p = sets[i].params;
            return 0;
        }
    }
    printf("unknown parameter set %s\n", argv[1]);
    return -1;
#else
    (void) argc;
    (void) argv;
    *p = NULL;
    return 0;
#endif
}

#endif
