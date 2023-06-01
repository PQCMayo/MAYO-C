// SPDX-License-Identifier: Apache-2.0

#ifndef MEM_H
#define MEM_H
#include <stddef.h>
#include <stdint.h>

#if defined(__GNUC__) || defined(__clang__)
#define BSWAP32(i) __builtin_bswap32((i))
#else
#define BSWAP32(i) ((((i) >> 24) & 0xff) | (((i) >> 8) & 0xff00) | (((i) & 0xff00) << 8) | ((i) << 24))
#endif

// a > b -> b - a is negative
// returns 0xFFFFFFFF if true, 0x00000000 if false
static inline uint32_t ct_is_greater_than(int a, int b) {
    int32_t diff = b - a;
    return (uint32_t) (diff >> (8*sizeof(uint32_t)-1));
}

// if a == b -> 0x00000000, else 0xFFFFFFFF
static inline uint32_t ct_compare_32(int a, int b) {
    return (uint32_t)((-(int32_t)(a ^ b)) >> (8*sizeof(uint32_t)-1));
}

// if a == b -> 0x00, else 0xFF
static inline unsigned char ct_compare_8(unsigned char a, unsigned char b) {
    return (int8_t)((-(int32_t)(a ^ b)) >> (8*sizeof(uint32_t)-1));
}

/**
 * Clears and frees allocated memory.
 * 
 * @param[out] mem Memory to be cleared and freed.
 * @param size Size of memory to be cleared and freed.
 */
void mayo_secure_free(void *mem, size_t size);

/**
 * Clears memory.
 * 
 * @param[out] mem Memory to be cleared.
 * @param size Size of memory to be cleared.
 */
void mayo_secure_clear(void *mem, size_t size);

#endif