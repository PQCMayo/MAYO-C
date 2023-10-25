// SPDX-License-Identifier: Apache-2.0

#ifndef BITSLICED_ARITHMETIC_96_H
#define BITSLICED_ARITHMETIC_96_H

#include <stdint.h>
#include <mayo.h>
#include <immintrin.h>

// This implements arithmetic for bitsliced vectors of 128 field elements in Z_2[x]/(x^4+x+1)

static
inline void bitsliced_96_vec_copy(const uint32_t *in, uint32_t *out) {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
    out[3] = in[3];
    out[4] = in[4];
    out[5] = in[5];
    out[6] = in[6];
    out[7] = in[7];
    out[8] = in[8];
    out[9] = in[9];
    out[10] = in[10];
    out[11] = in[11];
}

static
inline void bitsliced_96_vec_add(const uint32_t *_in, uint32_t *_acc) {
    const __m128i* in = (const __m128i*) _in;
    __m128i* acc = (__m128i*) _acc;
    acc[0] ^= in[0];
    acc[1] ^= in[1];
    acc[2] ^= in[2];
}

static 
inline void bitsliced_96_vec_mul_add_x(const uint32_t *__in, uint32_t *__acc) {
    const __m128i* in = (const __m128i*) __in;
    __m128i* acc = (__m128i*) __acc;
    
    __m128i in09101100 = _mm_alignr_epi8(in[0], in[2], 1*4);
    __m128i in01020304 = _mm_alignr_epi8(in[1], in[0], 1*4);
    __m128i in05060708 = _mm_alignr_epi8(in[2], in[1], 1*4);
    __m128i inXXXXXX09 = _mm_slli_si128(in09101100, 3*4);
    __m128i in1011XXXX = _mm_srli_si128(in[2], 2*4);

    acc[0] ^= in09101100 ^ inXXXXXX09;
    acc[1] ^= in01020304 ^ in1011XXXX;
    acc[2] ^= in05060708;
}

static 
inline void bitsliced_96_vec_mul_add_x_inv(const uint32_t *_in, uint32_t *_acc) {
    const __m128i* in = (const __m128i*) _in;
    __m128i* acc = (__m128i*) _acc;

    const __m128i zero_mask = _mm_setr_epi32(-1, -1, -1, 0);
    __m128i in03040506 = _mm_alignr_epi8(in[1], in[0], 3*4);
    __m128i in000102XX = _mm_and_si128(in[0], zero_mask);
    __m128i in07080910 = _mm_alignr_epi8(in[2], in[1], 3*4);
    __m128i in11000102 = _mm_alignr_epi8(in[0], in[2], 3*4);

    acc[0] ^= in000102XX ^ in03040506;
    acc[1] ^= in07080910;
    acc[2] ^= in11000102;
}

static 
inline void bitsliced_96_vec_mul_add(const uint32_t *_in, unsigned char a, uint32_t *_acc) {
    const __m128i* in = (const __m128i*) _in;
    __m128i* acc = (__m128i*) _acc;

    const __m128i lut_a0 = _mm_setr_epi8(0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1);
    const __m128i lut_a1 = _mm_setr_epi8(0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1);
    const __m128i lut_a2 = _mm_setr_epi8(0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1);
    const __m128i lut_a3 = _mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1);

    __m128i aa = _mm_set1_epi8(a);
    __m128i a0 = _mm_shuffle_epi8(lut_a0, aa);
    __m128i a1 = _mm_shuffle_epi8(lut_a1, aa);
    __m128i a2 = _mm_shuffle_epi8(lut_a2, aa);
    __m128i a3 = _mm_shuffle_epi8(lut_a3, aa);

    __m128i inrot[3];
    inrot[0] = _mm_alignr_epi8(in[1], in[0], 3*4);
    inrot[1] = _mm_alignr_epi8(in[2], in[1], 3*4);
    inrot[2] = _mm_alignr_epi8(in[0], in[2], 3*4);

    __m128i inrot2[3];
    inrot2[0] = _mm_alignr_epi8(in[2], in[1], 2*4);
    inrot2[1] = _mm_alignr_epi8(in[0], in[2], 2*4);
    inrot2[2] = _mm_alignr_epi8(in[1], in[0], 2*4);

    __m128i inrot3[3];
    inrot3[0] = _mm_alignr_epi8(in[0], in[2], 1*4);
    inrot3[1] = _mm_alignr_epi8(in[1], in[0], 1*4);
    inrot3[2] = _mm_alignr_epi8(in[2], in[1], 1*4);

    acc[0] ^= a0 & in[0];
    acc[1] ^= a0 & in[1];
    acc[2] ^= a0 & in[2];

    const __m128i andmask1 = _mm_setr_epi32(0, 0, 0, -1);
    const __m128i andmask2 = _mm_setr_epi32(-1, -1, 0, 0);
    const __m128i andmask3 = _mm_setr_epi32(-1, 0, 0, 0);
    acc[0] ^= a1 & (inrot3[0] ^ (inrot2[0] & andmask1));
    acc[1] ^= a1 & (inrot3[1] ^ (inrot2[1] & andmask2));
    acc[2] ^= a1 & (inrot3[2]);

    acc[0] ^= a2 & (inrot2[0] ^ (inrot[0] & andmask1));
    acc[1] ^= a2 & (inrot2[1] ^ (inrot[1]));
    acc[2] ^= a2 & (inrot2[2] ^ (inrot[2] & andmask3));

    acc[0] ^= a3 & (inrot[0] ^ (in[0] & andmask1));
    acc[1] ^= a3 & (inrot[1] ^ (in[1]));
    acc[2] ^= a3 & (inrot[2] ^ (in[2]));
}

static 
inline void bitsliced_96_multiply_bins(uint32_t *bins, uint32_t *out) {
    bitsliced_96_vec_mul_add_x_inv(bins +  5 * 12, bins +  10 * 12);
    bitsliced_96_vec_mul_add_x(bins + 11 * 12, bins + 12 * 12);
    bitsliced_96_vec_mul_add_x_inv(bins +  10 * 12, bins +  7 * 12);
    bitsliced_96_vec_mul_add_x(bins + 12 * 12, bins +  6 * 12);
    bitsliced_96_vec_mul_add_x_inv(bins +  7 * 12, bins +  14 * 12);
    bitsliced_96_vec_mul_add_x(bins +  6 * 12, bins +  3 * 12);
    bitsliced_96_vec_mul_add_x_inv(bins +  14 * 12, bins +  15 * 12);
    bitsliced_96_vec_mul_add_x(bins +  3 * 12, bins +  8 * 12);
    bitsliced_96_vec_mul_add_x_inv(bins +  15 * 12, bins +  13 * 12);
    bitsliced_96_vec_mul_add_x(bins +  8 * 12, bins +  4 * 12);
    bitsliced_96_vec_mul_add_x_inv(bins +  13 * 12, bins +  9 * 12);
    bitsliced_96_vec_mul_add_x(bins +  4 * 12, bins +  2 * 12);
    bitsliced_96_vec_mul_add_x_inv(bins +   9 * 12, bins +  1 * 12);
    bitsliced_96_vec_mul_add_x(bins +  2 * 12, bins +  1 * 12);
    bitsliced_96_vec_copy(bins + 12, out);
}

#endif

static
inline void mayo_3_P1_times_O(const uint32_t *_P1, const unsigned char *O, uint32_t *_acc) {
    const __m128i* P1 = (const __m128i*) _P1;
    __m128i* acc = (__m128i*) _acc;

    const __m128i lut_a0 = _mm_setr_epi8(0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1);
    const __m128i lut_a1 = _mm_setr_epi8(0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1);
    const __m128i lut_a2 = _mm_setr_epi8(0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1);
    const __m128i lut_a3 = _mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i andmask1 = _mm_setr_epi32(0, 0, 0, -1);
    const __m128i andmask2 = _mm_setr_epi32(-1, -1, 0, 0);
    const __m128i andmask3 = _mm_setr_epi32(-1, 0, 0, 0);

    #define MAYO_POS ( 3*(r*V_MAX + c - (r)*(r+1)/2) )

    for (int c = 0; c < V_MAX; c++) {
        for (int k = 0; k < O_MAX; k += 10) {

            #define PART(x) \
            __m128i aa##x = _mm_set1_epi8(O[c * O_MAX + k + x]); \
            __m128i a0##x = _mm_shuffle_epi8(lut_a0, aa##x); \
            __m128i a1##x = _mm_shuffle_epi8(lut_a1, aa##x); \
            __m128i a2##x = _mm_shuffle_epi8(lut_a2, aa##x); \
            __m128i a3##x = _mm_shuffle_epi8(lut_a3, aa##x); \

            PART(0)
            PART(1)
            PART(2)
            PART(3)
            PART(4)
            PART(5)
            PART(6)
            PART(7)
            PART(8)
            PART(9)

            #undef PART

            for (int r = 0; r <= c; r++) {

                int mp = MAYO_POS;

                __m128i in[3];
                in[0] = P1[mp + 0];
                in[1] = P1[mp + 1];
                in[2] = P1[mp + 2];

                __m128i inrot[3];
                inrot[0] = _mm_alignr_epi8(in[1], in[0], 3*4);
                inrot[1] = _mm_alignr_epi8(in[2], in[1], 3*4);
                inrot[2] = _mm_alignr_epi8(in[0], in[2], 3*4);

                __m128i inrot2[3];
                inrot2[0] = _mm_alignr_epi8(in[2], in[1], 2*4);
                inrot2[1] = _mm_alignr_epi8(in[0], in[2], 2*4);
                inrot2[2] = _mm_alignr_epi8(in[1], in[0], 2*4);

                __m128i inrot3[3];
                inrot3[0] = _mm_alignr_epi8(in[0], in[2], 1*4);
                inrot3[1] = _mm_alignr_epi8(in[1], in[0], 1*4);
                inrot3[2] = _mm_alignr_epi8(in[2], in[1], 1*4);

                __m128i X1 = inrot3[0] ^ (inrot2[0] & andmask1);
                __m128i X2 = inrot3[1] ^ (inrot2[1] & andmask2);
                __m128i X3 = inrot2[0] ^ (inrot[0] & andmask1);
                __m128i X4 = inrot2[1] ^ inrot[1];
                __m128i X5 = inrot2[2] ^ (inrot[2] & andmask3);
                __m128i X6 = inrot[0] ^ (in[0] & andmask1);
                __m128i X7 = inrot[1] ^ (in[1]);
                __m128i X8 = inrot[2] ^ (in[2]);

                #define PART(x) \
                acc[3 * (r * O_MAX + k + x) + 0] ^= a0##x & in[0]; \
                acc[3 * (r * O_MAX + k + x) + 1] ^= a0##x & in[1]; \
                acc[3 * (r * O_MAX + k + x) + 2] ^= a0##x & in[2]; \
                acc[3 * (r * O_MAX + k + x) + 0] ^= a1##x & X1; \
                acc[3 * (r * O_MAX + k + x) + 1] ^= a1##x & X2; \
                acc[3 * (r * O_MAX + k + x) + 2] ^= a1##x & (inrot3[2]); \
                acc[3 * (r * O_MAX + k + x) + 0] ^= a2##x & X3; \
                acc[3 * (r * O_MAX + k + x) + 1] ^= a2##x & X4; \
                acc[3 * (r * O_MAX + k + x) + 2] ^= a2##x & X5; \
                acc[3 * (r * O_MAX + k + x) + 0] ^= a3##x & X6; \
                acc[3 * (r * O_MAX + k + x) + 1] ^= a3##x & X7; \
                acc[3 * (r * O_MAX + k + x) + 2] ^= a3##x & X8; \

                PART(0)
                PART(1)
                PART(2)
                PART(3)
                PART(4)
                PART(5)
                PART(6)
                PART(7)
                PART(8)
                PART(9)
                #undef PART
            }
        }
    }

    #undef MAYO_POS
}

static
inline void mayo_3_P1P1t_times_O(const uint32_t *_P1P1t, const unsigned char *O, uint32_t *_acc) {
    const __m128i* P1P1t = (const __m128i*) _P1P1t;
    __m128i* acc = (__m128i*) _acc;

    const __m128i lut_a0 = _mm_setr_epi8(0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1);
    const __m128i lut_a1 = _mm_setr_epi8(0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1);
    const __m128i lut_a2 = _mm_setr_epi8(0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1);
    const __m128i lut_a3 = _mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i andmask1 = _mm_setr_epi32(0, 0, 0, -1);
    const __m128i andmask2 = _mm_setr_epi32(-1, -1, 0, 0);
    const __m128i andmask3 = _mm_setr_epi32(-1, 0, 0, 0);

    #define MAYO_POS ( 3*(r*V_MAX + c) )

    for (int c = 0; c < V_MAX; c++) {
        for (int k = 0; k < O_MAX; k += 10) {

            #define PART(x) \
            __m128i aa##x = _mm_set1_epi8(O[c * O_MAX + k + x]); \
            __m128i a0##x = _mm_shuffle_epi8(lut_a0, aa##x); \
            __m128i a1##x = _mm_shuffle_epi8(lut_a1, aa##x); \
            __m128i a2##x = _mm_shuffle_epi8(lut_a2, aa##x); \
            __m128i a3##x = _mm_shuffle_epi8(lut_a3, aa##x); \

            PART(0)
            PART(1)
            PART(2)
            PART(3)
            PART(4)
            PART(5)
            PART(6)
            PART(7)
            PART(8)
            PART(9)

            #undef PART

            for (int r = 0; r < V_MAX; r++) {

                int mp = MAYO_POS;

                __m128i in[3];
                in[0] = P1P1t[mp + 0];
                in[1] = P1P1t[mp + 1];
                in[2] = P1P1t[mp + 2];

                __m128i inrot[3];
                inrot[0] = _mm_alignr_epi8(in[1], in[0], 3*4);
                inrot[1] = _mm_alignr_epi8(in[2], in[1], 3*4);
                inrot[2] = _mm_alignr_epi8(in[0], in[2], 3*4);

                __m128i inrot2[3];
                inrot2[0] = _mm_alignr_epi8(in[2], in[1], 2*4);
                inrot2[1] = _mm_alignr_epi8(in[0], in[2], 2*4);
                inrot2[2] = _mm_alignr_epi8(in[1], in[0], 2*4);

                __m128i inrot3[3];
                inrot3[0] = _mm_alignr_epi8(in[0], in[2], 1*4);
                inrot3[1] = _mm_alignr_epi8(in[1], in[0], 1*4);
                inrot3[2] = _mm_alignr_epi8(in[2], in[1], 1*4);

                __m128i X1 = inrot3[0] ^ (inrot2[0] & andmask1);
                __m128i X2 = inrot3[1] ^ (inrot2[1] & andmask2);
                __m128i X3 = inrot2[0] ^ (inrot[0] & andmask1);
                __m128i X4 = inrot2[1] ^ inrot[1];
                __m128i X5 = inrot2[2] ^ (inrot[2] & andmask3);
                __m128i X6 = inrot[0] ^ (in[0] & andmask1);
                __m128i X7 = inrot[1] ^ (in[1]);
                __m128i X8 = inrot[2] ^ (in[2]);

                #define PART(x) \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 0], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 0]) ^ (a0##x & in[0])); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 1], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 1]) ^ (a0##x & in[1])); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 2], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 2]) ^ (a0##x & in[2])); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 0], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 0]) ^ (a1##x & X1)); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 1], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 1]) ^ (a1##x & X2)); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 2], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 2]) ^ (a1##x & (inrot3[2]))); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 0], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 0]) ^ (a2##x & X3)); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 1], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 1]) ^ (a2##x & X4)); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 2], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 2]) ^ (a2##x & X5)); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 0], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 0]) ^ (a3##x & X6)); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 1], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 1]) ^ (a3##x & X7)); \
                _mm_storeu_si128(&acc[3 * (r * O_MAX + k + x) + 2], _mm_loadu_si128(&acc[3 * (r * O_MAX + k + x) + 2]) ^ (a3##x & X8)); \

                PART(0)
                PART(1)
                PART(2)
                PART(3)
                PART(4)
                PART(5)
                PART(6)
                PART(7)
                PART(8)
                PART(9)
                #undef PART
            }
        }
    }

    #undef MAYO_POS
}


static
inline void mayo_3_P1_times_Vt(const uint32_t *_P1, const unsigned char *V, uint32_t *_acc) {
    const __m128i* P1 = (const __m128i*) _P1;
    __m128i* acc = (__m128i*) _acc;

    const __m128i lut_a0 = _mm_setr_epi8(0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1);
    const __m128i lut_a1 = _mm_setr_epi8(0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1);
    const __m128i lut_a2 = _mm_setr_epi8(0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1);
    const __m128i lut_a3 = _mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i andmask1 = _mm_setr_epi32(0, 0, 0, -1);
    const __m128i andmask2 = _mm_setr_epi32(-1, -1, 0, 0);
    const __m128i andmask3 = _mm_setr_epi32(-1, 0, 0, 0);

    #define MAYO_POS ( 3*(r*V_MAX + c - (r)*(r+1)/2) )

    for (int c = 0; c < V_MAX; c++) {
        for (int k = 0; k < K_MAX; k += 11) {

            #define PART(x) \
            __m128i aa##x = _mm_set1_epi8(V[(k + x) * V_MAX + c]); \
            __m128i a0##x = _mm_shuffle_epi8(lut_a0, aa##x); \
            __m128i a1##x = _mm_shuffle_epi8(lut_a1, aa##x); \
            __m128i a2##x = _mm_shuffle_epi8(lut_a2, aa##x); \
            __m128i a3##x = _mm_shuffle_epi8(lut_a3, aa##x); \

            PART(0)
            PART(1)
            PART(2)
            PART(3)
            PART(4)
            PART(5)
            PART(6)
            PART(7)
            PART(8)
            PART(9)
            PART(10)

            #undef PART

            for (int r = 0; r <= c; r++) {

                int mp = MAYO_POS;

                __m128i in[3];
                in[0] = P1[mp + 0];
                in[1] = P1[mp + 1];
                in[2] = P1[mp + 2];

                __m128i inrot[3];
                inrot[0] = _mm_alignr_epi8(in[1], in[0], 3*4);
                inrot[1] = _mm_alignr_epi8(in[2], in[1], 3*4);
                inrot[2] = _mm_alignr_epi8(in[0], in[2], 3*4);

                __m128i inrot2[3];
                inrot2[0] = _mm_alignr_epi8(in[2], in[1], 2*4);
                inrot2[1] = _mm_alignr_epi8(in[0], in[2], 2*4);
                inrot2[2] = _mm_alignr_epi8(in[1], in[0], 2*4);

                __m128i inrot3[3];
                inrot3[0] = _mm_alignr_epi8(in[0], in[2], 1*4);
                inrot3[1] = _mm_alignr_epi8(in[1], in[0], 1*4);
                inrot3[2] = _mm_alignr_epi8(in[2], in[1], 1*4);

                __m128i X1 = inrot3[0] ^ (inrot2[0] & andmask1);
                __m128i X2 = inrot3[1] ^ (inrot2[1] & andmask2);
                __m128i X3 = inrot2[0] ^ (inrot[0] & andmask1);
                __m128i X4 = inrot2[1] ^ inrot[1];
                __m128i X5 = inrot2[2] ^ (inrot[2] & andmask3);
                __m128i X6 = inrot[0] ^ (in[0] & andmask1);
                __m128i X7 = inrot[1] ^ (in[1]);
                __m128i X8 = inrot[2] ^ (in[2]);

                #define PART(x) \
                acc[3 * (r * K_MAX + k + x) + 0] ^= a0##x & in[0]; \
                acc[3 * (r * K_MAX + k + x) + 1] ^= a0##x & in[1]; \
                acc[3 * (r * K_MAX + k + x) + 2] ^= a0##x & in[2]; \
                acc[3 * (r * K_MAX + k + x) + 0] ^= a1##x & X1; \
                acc[3 * (r * K_MAX + k + x) + 1] ^= a1##x & X2; \
                acc[3 * (r * K_MAX + k + x) + 2] ^= a1##x & (inrot3[2]); \
                acc[3 * (r * K_MAX + k + x) + 0] ^= a2##x & X3; \
                acc[3 * (r * K_MAX + k + x) + 1] ^= a2##x & X4; \
                acc[3 * (r * K_MAX + k + x) + 2] ^= a2##x & X5; \
                acc[3 * (r * K_MAX + k + x) + 0] ^= a3##x & X6; \
                acc[3 * (r * K_MAX + k + x) + 1] ^= a3##x & X7; \
                acc[3 * (r * K_MAX + k + x) + 2] ^= a3##x & X8; \

                PART(0)
                PART(1)
                PART(2)
                PART(3)
                PART(4)
                PART(5)
                PART(6)
                PART(7)
                PART(8)
                PART(9)
                PART(10)
                #undef PART
            }
        }
    }

    #undef MAYO_POS
}