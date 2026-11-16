// SPDX-License-Identifier: Apache-2.0

#ifndef M_VEC_ARITHMETIC_H
#define M_VEC_ARITHMETIC_H

#include <stdint.h>
#include <mayo.h>
#include <simple_arithmetic.h>

// This implements arithmetic for vectors of m field elements in Z_2[x]/(x^4+x+1),
// nibble-packed into m_vec_limbs 64-bit limbs.

/* The number of limbs the loops run over. With fixed parameters it is a compile-time
   constant, which lets the compiler unroll; the m_vec_limbs arguments then only
   matter with dynamic parameters. */
#ifdef ENABLE_PARAMS_DYNAMIC
#define MAYO_MVL(m_vec_limbs) (m_vec_limbs)
#else
#define MAYO_MVL(m_vec_limbs) M_VEC_LIMBS_MAX
#endif

static
inline void m_vec_copy (int m_vec_limbs, const uint64_t *in, uint64_t *out) {
    (void) m_vec_limbs;
    for (size_t i = 0; i < (size_t) MAYO_MVL(m_vec_limbs); i++)
    {
        out[i] = in[i];
    }
}

static
inline void m_vec_add (int m_vec_limbs, const uint64_t *in, uint64_t *acc) {
    (void) m_vec_limbs;
    for (size_t i = 0; i < (size_t) MAYO_MVL(m_vec_limbs); i++)
    {
        acc[i] ^= in[i];
    }
}

static
inline void m_vec_mul_add (int m_vec_limbs, const uint64_t *in, unsigned char a, uint64_t *acc) {
    (void) m_vec_limbs;
    uint32_t tab = mul_table(a);

    uint64_t lsb_ask = 0x1111111111111111ULL;

    for(int i=0; i < MAYO_MVL(m_vec_limbs) ;i++){
        acc[i] ^= ( in[i]       & lsb_ask) * (tab & 0xff)
                ^ ((in[i] >> 1) & lsb_ask) * ((tab >> 8)  & 0xf)
                ^ ((in[i] >> 2) & lsb_ask) * ((tab >> 16) & 0xf)
                ^ ((in[i] >> 3) & lsb_ask) * ((tab >> 24) & 0xf);
    }
}

inline
static void m_vec_mul_add_x (int m_vec_limbs, const uint64_t *in, uint64_t *acc) {
    (void) m_vec_limbs;
    uint64_t mask_msb = 0x8888888888888888ULL;
    for(int i=0; i < MAYO_MVL(m_vec_limbs); i++){
        uint64_t t = in[i] & mask_msb;
        acc[i] ^= ((in[i] ^ t) << 1) ^ ((t >> 3) * 3);
    }
}

inline
static void m_vec_mul_add_x_inv (int m_vec_limbs, const uint64_t *in, uint64_t *acc) {
    (void) m_vec_limbs;
    uint64_t mask_lsb = 0x1111111111111111ULL;
    for(int i=0; i < MAYO_MVL(m_vec_limbs); i++){
        uint64_t t = in[i] & mask_lsb;
        acc[i] ^= ((in[i] ^ t) >> 1) ^ (t * 9);
    }
}

static
inline void m_vec_multiply_bins (int m_vec_limbs, uint64_t *bins, uint64_t *out) {
    const int mvl = MAYO_MVL(m_vec_limbs);
    m_vec_mul_add_x_inv (m_vec_limbs, bins +  5 * mvl, bins + 10 * mvl);
    m_vec_mul_add_x     (m_vec_limbs, bins + 11 * mvl, bins + 12 * mvl);
    m_vec_mul_add_x_inv (m_vec_limbs, bins + 10 * mvl, bins +  7 * mvl);
    m_vec_mul_add_x     (m_vec_limbs, bins + 12 * mvl, bins +  6 * mvl);
    m_vec_mul_add_x_inv (m_vec_limbs, bins +  7 * mvl, bins + 14 * mvl);
    m_vec_mul_add_x     (m_vec_limbs, bins +  6 * mvl, bins +  3 * mvl);
    m_vec_mul_add_x_inv (m_vec_limbs, bins + 14 * mvl, bins + 15 * mvl);
    m_vec_mul_add_x     (m_vec_limbs, bins +  3 * mvl, bins +  8 * mvl);
    m_vec_mul_add_x_inv (m_vec_limbs, bins + 15 * mvl, bins + 13 * mvl);
    m_vec_mul_add_x     (m_vec_limbs, bins +  8 * mvl, bins +  4 * mvl);
    m_vec_mul_add_x_inv (m_vec_limbs, bins + 13 * mvl, bins +  9 * mvl);
    m_vec_mul_add_x     (m_vec_limbs, bins +  4 * mvl, bins +  2 * mvl);
    m_vec_mul_add_x_inv (m_vec_limbs, bins +  9 * mvl, bins +  1 * mvl);
    m_vec_mul_add_x     (m_vec_limbs, bins +  2 * mvl, bins +  1 * mvl);
    m_vec_copy          (m_vec_limbs, bins + mvl, out);
}

#endif
