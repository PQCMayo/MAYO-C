// SPDX-License-Identifier: Apache-2.0

#ifndef ARITHMETIC_COMMON_H
#define ARITHMETIC_COMMON_H

#include <mayo.h>
#include <arithmetic_64.h>
#include <arithmetic_96.h>
#include <arithmetic_128.h>
#include <stdalign.h>
#include <stdint.h>

#include <arm_neon.h>

const unsigned char __0_f[16] __attribute__((aligned(16))) = {
0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07, 0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f
};

const unsigned char __gf16_reduce[16] __attribute__((aligned(16))) = {
0x00,0x13,0x26,0x35,0x4c,0x5f,0x6a,0x79, 0x8b,0x98,0xad,0xbe,0xc7,0xd4,0xe1,0xf2
};

static inline
uint8x16_t _gf16v_mul_unpack_neon( uint8x16_t a0 , uint8x16_t b0 , uint8x16_t tab_reduce )
{
    uint8x16_t ab = vreinterpretq_u8_p8(vmulq_p8( a0 , b0 ));
    return ab^vqtbl1q_u8( tab_reduce , vshrq_n_u8(ab,4) );
}

static inline
uint8x16_t _gf16v_get_multab_neon( uint8x16_t b , uint8x16_t tab_reduce , uint8x16_t tab_0_f ) { return _gf16v_mul_unpack_neon(b,tab_0_f,tab_reduce); }

static inline
uint8x16_t gf16v_get_multab_neon( uint8_t b )
{
    uint8x16_t tab_reduce = vld1q_u8(__gf16_reduce);
    uint8x16_t tab_0_f = vld1q_u8(__0_f);

    uint8x16_t bb = vdupq_n_u8(b);
    return _gf16v_get_multab_neon(bb,tab_reduce,tab_0_f);
}

static inline
uint8x16_t _gf16_tbl_x2( uint8x16_t a , uint8x16_t tbl , uint8x16_t tblhi , uint8x16_t mask_f ) {
    // return vsliq_n_u8( vqtbl1q_u8( tbl , a&mask_f ) , vqtbl1q_u8( tbl , vshrq_n_u8( a , 4 ) ), 4 );
    return vqtbl1q_u8( tbl , a&mask_f ) ^ vqtbl1q_u8( tblhi , vshrq_n_u8( a , 4 ) );
}


#ifndef MAYO_VARIANT
static void m_multiply_bins(const int m_legs, uint64_t *bins, uint64_t *out) {

    m_vec_add(m_legs, bins + 15 * m_legs * 2, bins + 12 * m_legs * 2);
    m_vec_add(m_legs, bins + 15 * m_legs * 2, bins +  3 * m_legs * 2);

    m_vec_add(m_legs, bins + 14 * m_legs * 2, bins +  8 * m_legs * 2);
    m_vec_add(m_legs, bins + 14 * m_legs * 2, bins +  6 * m_legs * 2);

    m_vec_add(m_legs, bins + 13 * m_legs * 2, bins + 10 * m_legs * 2);
    m_vec_add(m_legs, bins + 13 * m_legs * 2, bins +  7 * m_legs * 2);

    m_vec_add(m_legs, bins + 12 * m_legs * 2, bins +  8 * m_legs * 2);
    m_vec_add(m_legs, bins + 12 * m_legs * 2, bins +  4 * m_legs * 2);

    m_vec_add(m_legs, bins + 11 * m_legs * 2, bins +  9 * m_legs * 2);
    m_vec_add(m_legs, bins + 11 * m_legs * 2, bins +  2 * m_legs * 2);

    m_vec_add(m_legs, bins + 10 * m_legs * 2, bins +  8 * m_legs * 2);
    m_vec_add(m_legs, bins + 10 * m_legs * 2, bins +  2 * m_legs * 2);

    m_vec_add(m_legs, bins + 9 * m_legs * 2, bins +  8 * m_legs * 2);
    m_vec_add(m_legs, bins + 9 * m_legs * 2, bins +  1 * m_legs * 2);

    m_vec_add(m_legs, bins + 7 * m_legs * 2, bins +  4 * m_legs * 2);
    m_vec_add(m_legs, bins + 7 * m_legs * 2, bins +  3 * m_legs * 2);

    m_vec_add(m_legs, bins + 6 * m_legs * 2, bins +  4 * m_legs * 2);
    m_vec_add(m_legs, bins + 6 * m_legs * 2, bins +  2 * m_legs * 2);

    m_vec_add(m_legs, bins + 5 * m_legs * 2, bins +  4 * m_legs * 2);
    m_vec_add(m_legs, bins + 5 * m_legs * 2, bins +  1 * m_legs * 2);

    m_vec_add(m_legs, bins + 3 * m_legs * 2, bins +  2 * m_legs * 2);
    m_vec_add(m_legs, bins + 3 * m_legs * 2, bins +  1 * m_legs * 2);

    m_vec_mul_add_x(m_legs, bins + 8 * m_legs * 2, bins + 4 * m_legs * 2);
    m_vec_mul_add_x(m_legs, bins + 4 * m_legs * 2, bins + 2 * m_legs * 2);
    m_vec_mul_add_x(m_legs, bins + 2 * m_legs * 2, bins + 1 * m_legs * 2);

    m_vec_copy(m_legs, bins + 1 * m_legs * 2, out);
}
#endif

// compute P * S^t = [ P1  P2 ] * [S1] = [P1*S1 + P2*S2]
//                   [  0  P3 ]   [S2]   [        P3*S2]
static inline void mayo_generic_m_calculate_PS(const uint64_t *P1, const uint64_t *P2, const uint64_t *P3, const unsigned char *S,
                              const int m, const int v, const int o, const int k, uint64_t *PS) {

    const int n = o + v;
#if defined(MAYO_VARIANT) && ((M_MAX == 64) || (M_MAX == 96) || M_MAX == 128)
    (void)m;
#else
    const int m_legs = m / 32;
#endif

    /* Old approach which is constant time but doesn't have to be
    unsigned char S1[V_MAX*K_MAX];
    unsigned char S2[O_MAX*K_MAX];
    unsigned char *s1_write = S1;
    unsigned char *s2_write = S2;
    for (int r=0; r < k; r++)
    {
        for (int c = 0; c < n; c++)
        {
            if(c < v){
                *(s1_write++) = S[r*n + c];
            } else {
                *(s2_write++) = S[r*n + c];
            }
        }
    }

    mul_add_m_upper_triangular_mat_x_mat_trans(m_legs, P1, S1, PS, v, v, k, 1); // P1 * S1
    mul_add_m_upper_triangular_mat_x_mat_trans(m_legs, P2, S2, PS, v, o, k, 0); // P2 * S2
    mul_add_m_upper_triangular_mat_x_mat_trans(m_legs, P3, S2, PS + v*k*m_legs*4, o, o, k, 1); // P3 * S2.
    */

    // use more stack efficient version for MAYO_3 and MAYO_5
    #if defined(PQM4) && N_MAX > 78
    uint64_t accumulator[M_MAX * N_MAX] = {0};
    int P1_used;
    int P3_used;
    for (int col = 0; col < k; col++) {
        for(unsigned int i = 0; i < sizeof(accumulator)/8; i++) {
            accumulator[i] = 0;
        }
        P1_used = 0;
        for (int row = 0; row < v; row++) {
            for (int j = row; j < v; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
                vec_add_64(P1 + P1_used * 4, accumulator + ( row * 16 + S[col * n + j] ) * 4 );
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
                vec_add_96(P1 + P1_used * 6, accumulator + ( row * 16 + S[col * n + j] ) * 6);
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
                vec_add_128(P1 + P1_used * 8, accumulator + ( row * 16 + S[col * n + j] ) * 8);
#else
                bitsliced_m_vec_add(m_legs, P1 + P1_used * m_legs * 2, accumulator + ( row * 16 + S[col * n + j] )*m_legs * 2 );
#endif
                P1_used ++;
            }

            for (int j = 0; j < o; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
                vec_add_64(P2 + (row * o + j)*4, accumulator + ( row * 16 + S[(col * n) + j + v] )*4 );
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
                vec_add_96(P2 + (row * o + j)*6, accumulator + ( row * 16 + S[(col * n) + j + v] )*6);
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
                vec_add_128(P2 + (row * o + j)*8, accumulator + ( row * 16 + S[(col * n) + j + v] )*8 );
#else
                bitsliced_m_vec_add(m_legs, P2 + (row * o + j)*m_legs * 2, accumulator + ( row * 16 + S[(col * n) + j + v] )*m_legs * 2 );
#endif
            }
        }

        P3_used = 0;
        for (int row = v; row < n; row++) {
            for (int j = row; j < n; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
                vec_add_64(P3 + P3_used * 4, accumulator + ( row * 16 + S[col * n + j] )*4);
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
                vec_add_96(P3 + P3_used * 6, accumulator + ( row * 16 + S[col * n + j] )*6);
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
                vec_add_128(P3 + P3_used * 8, accumulator + ( row * 16 + S[col * n + j] )*8);
#else
                bitsliced_m_vec_add(m_legs, P3 + P3_used * m_legs * 2, accumulator + ( row * 16 + S[col * n + j] )*m_legs * 2 );
#endif
                P3_used ++;
            }
        }

        for (int row = 0; row < n; row++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
           multiply_bins_64(accumulator + row * 16 * 4, PS + (row * k + col) * 4);
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
           multiply_bins_96(accumulator + row * 16 * 6, PS + (row * k + col) * 6);
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
           multiply_bins_128(accumulator + row * 16 * 8, PS + (row * k + col) * 8);
#else
           bitsliced_m_multiply_bins(m_legs, accumulator + row * 16 * m_legs * 2, PS + (row * k + col) * m_legs * 2);
#endif
        }
    }

    #else

    alignas (32) uint64_t accumulator[M_MAX * K_MAX * N_MAX] = {0};
    int P1_used = 0;
    for (int row = 0; row < v; row++) {
        for (int j = row; j < v; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 9)
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 4 );
#elif defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 4)
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 4 );
            vec_add_64(P1 + P1_used * 4, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 4 );
#elif defined(MAYO_VARIANT) && (M_MAX == 96) && (K_MAX == 11)
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 9) * 16 + S[9 * n + j] ) * 6 );
            vec_add_96(P1 + P1_used * 6, accumulator + ( (row * k + 10) * 16 + S[10 * n + j] ) * 6 );
#elif defined(MAYO_VARIANT) && (M_MAX == 128) && (K_MAX == 12)
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 9) * 16 + S[9 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 10) * 16 + S[10 * n + j] ) * 8 );
            vec_add_128(P1 + P1_used * 8, accumulator + ( (row * k + 11) * 16 + S[11 * n + j] ) * 8 );
#else
            for (int col = 0; col < k; col++) {
                m_vec_add(m_legs, P1 + P1_used * m_legs * 2, accumulator + ( (row * k + col) * 16 + S[col * n + j] )*m_legs * 2 );
            }
#endif
            P1_used ++;
        }


        for (int j = 0; j < o; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 9)
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 0) * 16 + S[(0 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 1) * 16 + S[(1 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 2) * 16 + S[(2 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 3) * 16 + S[(3 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 4) * 16 + S[(4 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 5) * 16 + S[(5 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 6) * 16 + S[(6 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 7) * 16 + S[(7 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 8) * 16 + S[(8 * n) + j + v] ) * 4 );
#elif defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 4)
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 0) * 16 + S[(0 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 1) * 16 + S[(1 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 2) * 16 + S[(2 * n) + j + v] ) * 4 );
            vec_add_64(P2 + (row * o + j) * 4, accumulator + ( (row * k + 3) * 16 + S[(3 * n) + j + v] ) * 4 );
#elif defined(MAYO_VARIANT) && (M_MAX == 96) && (K_MAX == 11)
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 0) * 16 + S[(0 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 1) * 16 + S[(1 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 2) * 16 + S[(2 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 3) * 16 + S[(3 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 4) * 16 + S[(4 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 5) * 16 + S[(5 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 6) * 16 + S[(6 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 7) * 16 + S[(7 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 8) * 16 + S[(8 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 9) * 16 + S[(9 * n) + j + v] ) * 6 );
            vec_add_96(P2 + (row * o + j) * 6, accumulator + ( (row * k + 10) * 16 + S[(10 * n) + j + v] ) * 6 );
#elif defined(MAYO_VARIANT) && (M_MAX == 128) && (K_MAX == 12)
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 0) * 16 + S[(0 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 1) * 16 + S[(1 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 2) * 16 + S[(2 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 3) * 16 + S[(3 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 4) * 16 + S[(4 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 5) * 16 + S[(5 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 6) * 16 + S[(6 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 7) * 16 + S[(7 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 8) * 16 + S[(8 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 9) * 16 + S[(9 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 10) * 16 + S[(10 * n) + j + v] ) * 8 );
            vec_add_128(P2 + (row * o + j) * 8, accumulator + ( (row * k + 11) * 16 + S[(11 * n) + j + v] ) * 8 );
#else
            for (int col = 0; col < k; col++) {
                m_vec_add(m_legs, P2 + (row * o + j)*m_legs * 2, accumulator + ( (row * k + col) * 16 + S[(col * n) + j + v] )*m_legs * 2 );
            }
#endif
        }
    }

    int P3_used = 0;
    for (int row = v; row < n; row++) {
        for (int j = row; j < n; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 9)
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 4 );
#elif defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 4)
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 4 );
            vec_add_64(P3 + P3_used * 4, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 4 );
#elif defined(MAYO_VARIANT) && (M_MAX == 96) && (K_MAX == 11)
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 9) * 16 + S[9 * n + j] ) * 6 );
            vec_add_96(P3 + P3_used * 6, accumulator + ( (row * k + 10) * 16 + S[10 * n + j] ) * 6 );
#elif defined(MAYO_VARIANT) && (M_MAX == 128) && (K_MAX == 12)
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 9) * 16 + S[9 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 10) * 16 + S[10 * n + j] ) * 8 );
            vec_add_128(P3 + P3_used * 8, accumulator + ( (row * k + 11) * 16 + S[11 * n + j] ) * 8 );
#else
            for (int col = 0; col < k; col++) {
                m_vec_add(m_legs, P3 + P3_used * m_legs * 2, accumulator + ( (row * k + col) * 16 + S[col * n + j] )*m_legs * 2 );
            }
#endif
            P3_used ++;
        }
    }

    // multiply stuff according to the bins of the accumulator and add to PS.
    int i = 0;
    while (i < n * k) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
        multiply_bins_64(accumulator + i * 16 * 4, PS + i * 4);
        i++;
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
        multiply_bins_96(accumulator + i * 16 * 6, PS + i * 6);
        i++;
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
        multiply_bins_128(accumulator + i * 16 * 8, PS + i * 8);
        i++;
#else
        m_multiply_bins(m/32, accumulator + i * 16 * (m/32) * 2, PS + i * (m/32) * 2);
        i++;
#endif
    }

    #endif
}


static inline void mayo_generic_m_calculate_SPS(const uint64_t *PS, const unsigned char *S, int m, int k, int  n, uint64_t *SPS){
    alignas (32) uint64_t accumulator[M_MAX*K_MAX*K_MAX] = {0};
    #if !defined(MAYO_VARIANT)
    const int m_legs = m/32;
    #else
    (void) m;
    #endif
    for (int row = 0; row < k; row++) {
        for (int j = 0; j < n; j++) {
            for (int col = 0; col < k; col += 1) {
                #if defined(MAYO_VARIANT) && (M_MAX == 64)
                    vec_add_64(PS + (j * k + col) * 4, accumulator + ( (row * k + col) * 16 + S[row * n + j] ) * 4 );
                #elif defined(MAYO_VARIANT) && (M_MAX == 96)
                    vec_add_96(PS + (j * k + col) * 6, accumulator + ( (row * k + col) * 16 + S[row * n + j] ) * 6 );
                #elif defined(MAYO_VARIANT) && (M_MAX == 128)
                    vec_add_128(PS + (j * k + col) * 8, accumulator + ( (row * k + col) * 16 + S[row * n + j] ) * 8 );
                #else
                    m_vec_add(m_legs, PS + (j * k + col) * m_legs * 2, accumulator + ( (row * k + col) * 16 + S[row * n + j] )*m_legs * 2 );
                #endif
            }
        }
    }

    // multiply stuff according to the bins of the accumulator and add to PS.
    int i = 0;
    while (i < k*k) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
        multiply_bins_64(accumulator + i * 16 * 4, SPS + i * 4);
        i++;
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
        multiply_bins_96(accumulator + i * 16 * 6, SPS + i * 6);
        i++;
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
        multiply_bins_128(accumulator + i * 16 * 8, SPS + i * 8);
        i++;
#else
        m_multiply_bins(m_legs, accumulator + i * 16 * m_legs * 2, SPS + i * m_legs * 2);
        i++;
#endif
    }
}


// multiplies m (possibly upper triangular) matrices with a single matrix and adds result to acc
static inline void mul_add_m_upper_triangular_mat_x_mat(int m_legs, const uint64_t *bs_mat, const unsigned char *mat, uint64_t *acc, int bs_mat_rows, int bs_mat_cols, int mat_cols, int triangular) {
    // assumes triangular == 1, i.e. bs_mat_cols == bs_mat_rows
    for (int k = 0; k < mat_cols; k++) {
        for (int c = 0; c < bs_mat_cols; c++) {
            uint8x16_t tbl = gf16v_get_multab_neon(mat[c * mat_cols + k]);
            uint8x16_t tblhi = tbl << 4;

            for (int r = 0; r <= c; r++) {
                int pos = r*(bs_mat_cols*2-r+1)/2 + (c - r);

                uint8_t *_acc = (uint8_t*)acc + (m_legs * 16) * (r * mat_cols + k);
                uint8_t *_bs = (uint8_t*)bs_mat + (m_legs * 16) * pos;
                for (int p = 0; p < m_legs * 16; p += 16) {
                    uint8x16_t t = vld1q_u8(_acc + p);
                    uint8x16_t a = vld1q_u8(_bs + p);

                    uint8x16_t mask_f = vdupq_n_u8( 0xf );
                    t ^= _gf16_tbl_x2(a, tbl, tblhi, mask_f);

                    vst1q_u8(_acc + p, t);
                }
            }
        }
    }
}

// multiplies m (possibly upper triangular) matrices with the transpose of a single matrix and adds result to acc
static inline void mul_add_m_upper_triangular_mat_x_mat_trans(int m_legs, const uint64_t *bs_mat, const unsigned char *mat, uint64_t *acc, int bs_mat_rows, int bs_mat_cols, int mat_rows, int triangular) {
    // assumes triangular == 1, i.e. bs_mat_cols == bs_mat_rows
    for (int k = 0; k < mat_rows; k++) {
        for (int c = 0; c < bs_mat_cols; c++) {
            uint8x16_t tbl = gf16v_get_multab_neon(mat[k * bs_mat_cols + c]);
            uint8x16_t tblhi = tbl << 4;

            for (int r = 0; r <= c; r++) {
                int pos = r*(bs_mat_cols*2-r+1)/2 + (c - r);

                uint8_t *_acc = (uint8_t*)acc + (m_legs * 16) * (r * mat_rows + k);
                uint8_t *_bs = (uint8_t*)bs_mat + (m_legs * 16) * pos;
                for (int p = 0; p < m_legs * 16; p += 16) {
                    uint8x16_t t = vld1q_u8(_acc + p);
                    uint8x16_t a = vld1q_u8(_bs + p);

                    uint8x16_t mask_f = vdupq_n_u8( 0xf );
                    t ^= _gf16_tbl_x2(a, tbl, tblhi, mask_f);

                    vst1q_u8(_acc + p, t);
                }
            }
        }
    }
}

// multiplies the transpose of a single matrix with m matrices and adds result to acc
static inline void mul_add_mat_trans_x_m_mat(int m_legs, const unsigned char *mat, const uint64_t *bs_mat, uint64_t *acc, int mat_rows, int mat_cols, int bs_mat_cols) {
    for (int r = 0; r < mat_cols; r++) {
        for (int c = 0; c < mat_rows; c++) {
            uint8x16_t tbl = gf16v_get_multab_neon(mat[c * mat_cols + r]);
            uint8x16_t tblhi = tbl << 4;

            int cols = bs_mat_cols * (m_legs * 32)/2;
            uint8_t *_acc = (uint8_t*)acc + r * cols;
            uint8_t *_bs = (uint8_t*)bs_mat + c * cols;

            for (int k = 0; k < cols; k += 16) {
                uint8x16_t t = vld1q_u8(_acc + k);
                uint8x16_t a = vld1q_u8(_bs + k);

                uint8x16_t mask_f = vdupq_n_u8( 0xf );
                t ^= _gf16_tbl_x2(a, tbl, tblhi, mask_f);

                vst1q_u8((uint8_t*)(acc) + (r * cols + k), t);
            }
        }
    }
}

// multiplies a single matrix with m matrices and adds result to acc
static inline void mul_add_mat_x_m_mat(int m_legs, const unsigned char *mat, const uint64_t *bs_mat, uint64_t *acc, int mat_rows, int mat_cols, int bs_mat_cols) {
    for (int r = 0; r < mat_rows; r++) {
        for (int c = 0; c < mat_cols; c++) {
            uint8x16_t tbl = gf16v_get_multab_neon(mat[r * mat_cols + c]);
            uint8x16_t tblhi = tbl << 4;

            int cols = bs_mat_cols * (m_legs * 32)/2;
            uint8_t *_acc = (uint8_t*)acc + r * cols;
            uint8_t *_bs = (uint8_t*)bs_mat + c * cols;

            for (int k = 0; k < cols; k += 16) {
                uint8x16_t t = vld1q_u8(_acc + k);
                uint8x16_t a = vld1q_u8(_bs + k);

                uint8x16_t mask_f = vdupq_n_u8( 0xf );
                t ^= _gf16_tbl_x2(a, tbl, tblhi, mask_f);

                vst1q_u8((uint8_t*)(acc) + (r * cols + k), t);
            }
        }
    }
}

#endif

