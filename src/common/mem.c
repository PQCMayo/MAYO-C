// SPDX-License-Identifier: Apache-2.0

#include <string.h>
#include <stdlib.h>
#include <stdint.h>

void mayo_secure_free(void *mem, size_t size) {
    if (mem) {
        typedef void *(*memset_t)(void *, int, size_t);
        static volatile memset_t memset_func = memset;
        memset_func(mem, 0, size);
        free(mem);
    }
}
void mayo_secure_clear(void *mem, size_t size) {
    typedef void *(*memset_t)(void *, int, size_t);
    static volatile memset_t memset_func = memset;
    memset_func(mem, 0, size);
}

volatile uint32_t uint32_t_blocker = 0;
volatile uint64_t uint64_t_blocker = 0;
volatile unsigned char unsigned_char_blocker = 0;
