
// SPDX-License-Identifier: Apache-2.0

#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include <stdint.h>
#include <mayo.h>
#include <stdint.h>
#include <stddef.h>

#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
#ifndef TARGET_BIG_ENDIAN
#define TARGET_BIG_ENDIAN
#endif
#endif

#if defined(MAYO_AVX) || defined(MAYO_NEON)
    #include <shuffle_arithmetic.h>
#elif defined(MAYO_M4)
    #include <m4_arithmetic.h>
#else
    #include <generic_arithmetic.h>
#endif

static
inline void vec_mul_add_u64(const int legs, const uint64_t *in, unsigned char a, uint64_t *acc) {
    uint32_t tab = mul_table(a);

    uint64_t lsb_ask = 0x1111111111111111ULL;

    for(int i=0; i < legs; i++){
        acc[i] ^= ( in[i]       & lsb_ask) * (tab & 0xff)
                ^ ((in[i] >> 1) & lsb_ask) * ((tab >> 8)  & 0xf)
                ^ ((in[i] >> 2) & lsb_ask) * ((tab >> 16) & 0xf)
                ^ ((in[i] >> 3) & lsb_ask) * ((tab >> 24) & 0xf);
    }
}

// Calculate Upper in KeyGen
#define m_upper MAYO_NAMESPACE(m_upper)
void m_upper(const mayo_params_t* p, const uint64_t *in, uint64_t *out, int size);

// Sample solution in Sign
#define sample_solution MAYO_NAMESPACE(sample_solution)
int sample_solution(const mayo_params_t *p, unsigned char *A, const unsigned char *y, const unsigned char *r, unsigned char *x, int k, int o, int m, int A_cols);

#endif

