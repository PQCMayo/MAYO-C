// SPDX-License-Identifier: Apache-2.0

#ifndef BITSLICED_ARITHMETIC_128_H
#define BITSLICED_ARITHMETIC_128_H

#include <stdint.h>
#include <mayo.h>
#include <immintrin.h>


// This implements arithmetic for bitsliced vectors of 128 field elements in Z_2[x]/(x^4+x+1)

static
inline void bitsliced_128_vec_copy(const uint64_t *in, uint64_t *out) {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
    out[3] = in[3];
    out[4] = in[4];
    out[5] = in[5];
    out[6] = in[6];
    out[7] = in[7];
}

static
inline void bitsliced_128_vec_sum(const uint64_t *_in1, const uint64_t *_in2, uint64_t *_out) {
    const __m256i* in1 = (const __m256i*) _in1;
    const __m256i* in2 = (const __m256i*) _in2;
    __m256i* out = (__m256i*) _out;
    out[0] = in1[0] ^ in2[0];
    out[1] = in1[1] ^ in2[1];
}

static
inline void bitsliced_128_vec_add(const uint64_t *_in, uint64_t *_acc) {
    const __m256i* in = (const __m256i*) _in;
    __m256i* acc = (__m256i*) _acc;
    acc[0] ^= in[0];
    acc[1] ^= in[1];
}

static 
inline void bitsliced_128_vec_mul_add_x(const uint64_t *_in, uint64_t *_acc) {
    const __m128i* in = (const __m128i*) _in;
    __m128i* acc = (__m128i*) _acc;
    acc[0] ^= in[3];
    acc[1] ^= in[0] ^ in[3];
    acc[2] ^= in[1];
    acc[3] ^= in[2];
}



static 
inline void bitsliced_128_vec_mul_add_x_inv(const uint64_t *_in, uint64_t *_acc) {
    const __m128i* in = (const __m128i*) _in;
    __m128i* acc = (__m128i*) _acc;
    acc[0] ^= in[0] ^ in[1];
    acc[1] ^= in[2];
    acc[2] ^= in[3];
    acc[3] ^= in[0];
}

static 
inline void bitsliced_128_vec_mul_add(const uint64_t *_in, unsigned char a, uint64_t *_acc) {
    const __m256i* in = (const __m256i*) _in;
    __m256i* acc = (__m256i*) _acc;

    const __m256i lut_a0 = _mm256_setr_epi8(
        0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1,
        0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1
    );
    const __m256i lut_a1 = _mm256_setr_epi8(
        0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1,
        0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, -1, -1
    );
    const __m256i lut_a2 = _mm256_setr_epi8(
        0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1,
        0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1
    );
    const __m256i lut_a3 = _mm256_setr_epi8(
        0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1,
        0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1
    );

    __m256i aa = _mm256_set1_epi8(a);
    __m256i a0 = _mm256_shuffle_epi8(lut_a0, aa);
    __m256i a1 = _mm256_shuffle_epi8(lut_a1, aa);
    __m256i a2 = _mm256_shuffle_epi8(lut_a2, aa);
    __m256i a3 = _mm256_shuffle_epi8(lut_a3, aa);

    __m256i xz = in[0] ^ _mm256_permute4x64_epi64(in[1], 0b01001110);
    __m256i yy = in[1] ^ _mm256_permute4x64_epi64(in[1], 0b01001110);

    // deg 0 term of a;
    acc[0] ^= a0 & in[0];
    acc[1] ^= a0 & in[1];
    
    // deg 1 term of a;
    acc[0] ^= a1 & _mm256_permute2x128_si256(xz,    in[1], 0x03);
    acc[1] ^= a1 & _mm256_permute2x128_si256(in[0], in[1], 0x21);

    // deg 2 term of a
    acc[0] ^= a2 & _mm256_blend_epi32(yy, in[1], 0x0F);
    acc[1] ^= a2 & _mm256_blend_epi32(xz, in[0], 0xF0);

    // deg 3 term of a
    acc[0] ^= a3 & _mm256_permute2x128_si256(xz, in[0], 0x13);
    acc[1] ^= a3 & _mm256_permute2x128_si256(xz, yy,    0x02);

    /*
     // deg 1 term of a;
    x = in[0] ^ in[3];
    acc[0] ^= a1 & in[3];
    acc[1] ^= a1 & x;
    acc[2] ^= a1 & in[1];
    acc[3] ^= a1 & in[2];

    // deg 2 term of a
    y = in[3] ^ in[2];
    acc[0] ^= a2 & in[2];
    acc[1] ^= a2 & y;
    acc[2] ^= a2 & x;
    acc[3] ^= a2 & in[1];

    // deg 3 term of a
    z = in[2] ^ in[1];
    acc[0] ^= a3 & in[1];
    acc[1] ^= a3 & z;
    acc[2] ^= a3 & y;
    acc[3] ^= a3 & x;
    */
}

static 
inline void bitsliced_128_multiply_bins(uint32_t *bins_32, uint32_t *out_32) {

    uint64_t *bins = (uint64_t *) bins_32;
    uint64_t *out = (uint64_t *) out_32;

    bitsliced_128_vec_mul_add_x_inv(bins +  5 * 8, bins +  10 * 8);
    bitsliced_128_vec_mul_add_x(bins + 11 * 8, bins + 12 * 8);
    bitsliced_128_vec_mul_add_x_inv(bins +  10 * 8, bins +  7 * 8);
    bitsliced_128_vec_mul_add_x(bins + 12 * 8, bins +  6 * 8);
    bitsliced_128_vec_mul_add_x_inv(bins +  7 * 8, bins +  14 * 8);
    bitsliced_128_vec_mul_add_x(bins +  6 * 8, bins +  3 * 8);
    bitsliced_128_vec_mul_add_x_inv(bins +  14 * 8, bins +  15 * 8);
    bitsliced_128_vec_mul_add_x(bins +  3 * 8, bins +  8 * 8);
    bitsliced_128_vec_mul_add_x_inv(bins +  15 * 8, bins +  13 * 8);
    bitsliced_128_vec_mul_add_x(bins +  8 * 8, bins +  4 * 8);
    bitsliced_128_vec_mul_add_x_inv(bins +  13 * 8, bins +  9 * 8);
    bitsliced_128_vec_mul_add_x(bins +  4 * 8, bins +  2 * 8);
    bitsliced_128_vec_mul_add_x_inv(bins +   9 * 8, bins +  1 * 8);
    bitsliced_128_vec_mul_add_x(bins +  2 * 8, bins +  1 * 8);
    bitsliced_128_vec_copy(bins + 8, out);
}

static
inline void mayo_5_P1_times_O(const uint32_t *_P1, const unsigned char *O, uint32_t *_acc) {

    const __m256i* P1 = (const __m256i*) _P1;
    __m256i* acc = (__m256i*) _acc;

    const __m256i mask11 = _mm256_set_epi64x(16,  16,  1 ,  1  );
    const __m256i mask12 = _mm256_set_epi64x(16,  16,  16,  16 );
    
    const __m256i mask21 = _mm256_set_epi64x(32,  32,  128, 128);
    const __m256i mask22 = _mm256_set_epi64x(32,  32,  8,   8);

    const __m256i mask31 = _mm256_set_epi64x( 4,  4,   64,  64 );
    const __m256i mask32 = _mm256_set_epi64x(64,  64,  64,  64 );

    const __m256i mask41 = _mm256_set_epi64x(8,   8,   32,  32 );
    const __m256i mask42 = _mm256_set_epi64x(128, 128, 32,  32 );

    const __m256i lut_a = _mm256_setr_epi8(
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1,
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1
    );

    #define MAYO_POS ( 2*(r*V_MAX + c - (r)*(r+1)/2) )

    for (int c = 0; c < V_MAX; c++) {
        int k;
        for (k = 0; k < (O_MAX/4)*4; k +=2) {

            __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k]));

            __m256i cmask1 = _mm256_cmpeq_epi64(aaaa & mask11, mask11);
            __m256i cmask2 = _mm256_cmpeq_epi64(aaaa & mask12, mask12);
            __m256i cmask3 = _mm256_cmpeq_epi64(aaaa & mask21, mask21);
            __m256i cmask4 = _mm256_cmpeq_epi64(aaaa & mask22, mask22);
            __m256i cmask5 = _mm256_cmpeq_epi64(aaaa & mask31, mask31);
            __m256i cmask6 = _mm256_cmpeq_epi64(aaaa & mask32, mask32);
            __m256i cmask7 = _mm256_cmpeq_epi64(aaaa & mask41, mask41);
            __m256i cmask8 = _mm256_cmpeq_epi64(aaaa & mask42, mask42);

            __m256i aaaa2 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 1]));

            __m256i cmask12 = _mm256_cmpeq_epi64(aaaa2 & mask11, mask11);
            __m256i cmask22 = _mm256_cmpeq_epi64(aaaa2 & mask12, mask12);
            __m256i cmask32 = _mm256_cmpeq_epi64(aaaa2 & mask21, mask21);
            __m256i cmask42 = _mm256_cmpeq_epi64(aaaa2 & mask22, mask22);
            __m256i cmask52 = _mm256_cmpeq_epi64(aaaa2 & mask31, mask31);
            __m256i cmask62 = _mm256_cmpeq_epi64(aaaa2 & mask32, mask32);
            __m256i cmask72 = _mm256_cmpeq_epi64(aaaa2 & mask41, mask41);
            __m256i cmask82 = _mm256_cmpeq_epi64(aaaa2 & mask42, mask42);

            for (int r = 0; r <= c; r++) {

                int mp = MAYO_POS;

                __m256i in0swap = _mm256_permute4x64_epi64(P1[mp+0], 0b01001110);
                __m256i in1swap = _mm256_permute4x64_epi64(P1[mp+1], 0b01001110);
                
                acc[2*(r * O_MAX + k) + 0] ^= (P1[mp+0] & cmask1) ^ (in0swap  & cmask3) 
                                            ^ (P1[mp+1] & cmask5) ^ (in1swap  & cmask7);
                acc[2*(r * O_MAX + k) + 1] ^= (P1[mp+1] & cmask2) ^ (in1swap  & cmask4)
                                            ^ (P1[mp+0] & cmask6) ^ (in0swap  & cmask8);

                acc[2*(r * O_MAX + k + 1) + 0] ^= (P1[mp+0] & cmask12) ^ (in0swap  & cmask32) 
                                                ^ (P1[mp+1] & cmask52) ^ (in1swap  & cmask72);
                acc[2*(r * O_MAX + k + 1) + 1] ^= (P1[mp+1] & cmask22) ^ (in1swap  & cmask42)
                                                ^ (P1[mp+0] & cmask62) ^ (in0swap  & cmask82);
            }
        }
    }

    #undef MAYO_POS
}


static
inline void mayo_5_P1P1t_times_O(const uint32_t *_P1P1t, const unsigned char *O, uint32_t *_acc) {

    const __m256i* P1P1t = (const __m256i*) _P1P1t;
    __m256i* acc = (__m256i*) _acc;

    const __m256i mask11 = _mm256_set_epi64x(16,  16,  1 ,  1  );
    const __m256i mask12 = _mm256_set_epi64x(16,  16,  16,  16 );
    
    const __m256i mask21 = _mm256_set_epi64x(32,  32,  128, 128);
    const __m256i mask22 = _mm256_set_epi64x(32,  32,  8,   8);

    const __m256i mask31 = _mm256_set_epi64x( 4,  4,   64,  64 );
    const __m256i mask32 = _mm256_set_epi64x(64,  64,  64,  64 );

    const __m256i mask41 = _mm256_set_epi64x(8,   8,   32,  32 );
    const __m256i mask42 = _mm256_set_epi64x(128, 128, 32,  32 );

    const __m256i lut_a = _mm256_setr_epi8(
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1,
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1
    );

    #define MAYO_POS ( 2*(r*V_MAX + c) )

    for (int c = 0; c < V_MAX; c++) {
        for (int k = 0; k < O_MAX; k +=2) {

            __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k]));

            __m256i cmask1 = _mm256_cmpeq_epi64(aaaa & mask11, mask11);
            __m256i cmask2 = _mm256_cmpeq_epi64(aaaa & mask12, mask12);
            __m256i cmask3 = _mm256_cmpeq_epi64(aaaa & mask21, mask21);
            __m256i cmask4 = _mm256_cmpeq_epi64(aaaa & mask22, mask22);
            __m256i cmask5 = _mm256_cmpeq_epi64(aaaa & mask31, mask31);
            __m256i cmask6 = _mm256_cmpeq_epi64(aaaa & mask32, mask32);
            __m256i cmask7 = _mm256_cmpeq_epi64(aaaa & mask41, mask41);
            __m256i cmask8 = _mm256_cmpeq_epi64(aaaa & mask42, mask42);

            __m256i aaaa2 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 1]));

            __m256i cmask12 = _mm256_cmpeq_epi64(aaaa2 & mask11, mask11);
            __m256i cmask22 = _mm256_cmpeq_epi64(aaaa2 & mask12, mask12);
            __m256i cmask32 = _mm256_cmpeq_epi64(aaaa2 & mask21, mask21);
            __m256i cmask42 = _mm256_cmpeq_epi64(aaaa2 & mask22, mask22);
            __m256i cmask52 = _mm256_cmpeq_epi64(aaaa2 & mask31, mask31);
            __m256i cmask62 = _mm256_cmpeq_epi64(aaaa2 & mask32, mask32);
            __m256i cmask72 = _mm256_cmpeq_epi64(aaaa2 & mask41, mask41);
            __m256i cmask82 = _mm256_cmpeq_epi64(aaaa2 & mask42, mask42);

            for (int r = 0; r < V_MAX; r++) {

                int mp = MAYO_POS;

                __m256i in0swap = _mm256_permute4x64_epi64(P1P1t[mp+0], 0b01001110);
                __m256i in1swap = _mm256_permute4x64_epi64(P1P1t[mp+1], 0b01001110);

                _mm256_storeu_si256(&acc[2*(r * O_MAX + k) + 0], _mm256_loadu_si256(&acc[2*(r * O_MAX + k) + 0]) ^ (
                    (P1P1t[mp+0] & cmask1) ^ (in0swap  & cmask3) 
                  ^ (P1P1t[mp+1] & cmask5) ^ (in1swap  & cmask7)
                ));
                _mm256_storeu_si256(&acc[2*(r * O_MAX + k) + 1], _mm256_loadu_si256(&acc[2*(r * O_MAX + k) + 1]) ^ (
                    (P1P1t[mp+1] & cmask2) ^ (in1swap  & cmask4) 
                  ^ (P1P1t[mp+0] & cmask6) ^ (in0swap  & cmask8)
                ));

                _mm256_storeu_si256(&acc[2*(r * O_MAX + k + 1) + 0], _mm256_loadu_si256(&acc[2*(r * O_MAX + k + 1) + 0]) ^ (
                    (P1P1t[mp+0] & cmask12) ^ (in0swap  & cmask32) 
                  ^ (P1P1t[mp+1] & cmask52) ^ (in1swap  & cmask72)
                ));
                _mm256_storeu_si256(&acc[2*(r * O_MAX + k + 1) + 1], _mm256_loadu_si256(&acc[2*(r * O_MAX + k + 1) + 1]) ^ (
                    (P1P1t[mp+1] & cmask22) ^ (in1swap  & cmask42) 
                  ^ (P1P1t[mp+0] & cmask62) ^ (in0swap  & cmask82)
                ));
            }
        }
    }

    #undef MAYO_POS
}


static
inline void mayo_5_P1_times_Vt(const uint32_t *_P1, const unsigned char *V, uint32_t *_acc) {

    const __m256i* P1 = (const __m256i*) _P1;
    __m256i* acc = (__m256i*) _acc;

    const __m256i mask11 = _mm256_set_epi64x(16,  16,  1 ,  1  );
    const __m256i mask12 = _mm256_set_epi64x(16,  16,  16,  16 );
    
    const __m256i mask21 = _mm256_set_epi64x(32,  32,  128, 128);
    const __m256i mask22 = _mm256_set_epi64x(32,  32,  8,   8);

    const __m256i mask31 = _mm256_set_epi64x( 4,  4,   64,  64 );
    const __m256i mask32 = _mm256_set_epi64x(64,  64,  64,  64 );

    const __m256i mask41 = _mm256_set_epi64x(8,   8,   32,  32 );
    const __m256i mask42 = _mm256_set_epi64x(128, 128, 32,  32 );

    const __m256i lut_a = _mm256_setr_epi8(
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1,
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1
    );

    #define MAYO_POS ( 2*(r*V_MAX + c - (r)*(r+1)/2) )

    for (int c = 0; c < V_MAX; c++) {
        for (int k = 0; k < K_MAX; k +=2) {

            __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[k * V_MAX + c]));

            __m256i cmask1 = _mm256_cmpeq_epi64(aaaa & mask11, mask11);
            __m256i cmask2 = _mm256_cmpeq_epi64(aaaa & mask12, mask12);
            __m256i cmask3 = _mm256_cmpeq_epi64(aaaa & mask21, mask21);
            __m256i cmask4 = _mm256_cmpeq_epi64(aaaa & mask22, mask22);
            __m256i cmask5 = _mm256_cmpeq_epi64(aaaa & mask31, mask31);
            __m256i cmask6 = _mm256_cmpeq_epi64(aaaa & mask32, mask32);
            __m256i cmask7 = _mm256_cmpeq_epi64(aaaa & mask41, mask41);
            __m256i cmask8 = _mm256_cmpeq_epi64(aaaa & mask42, mask42);

            __m256i aaaa2 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[(k+1) * V_MAX + c]));

            __m256i cmask12 = _mm256_cmpeq_epi64(aaaa2 & mask11, mask11);
            __m256i cmask22 = _mm256_cmpeq_epi64(aaaa2 & mask12, mask12);
            __m256i cmask32 = _mm256_cmpeq_epi64(aaaa2 & mask21, mask21);
            __m256i cmask42 = _mm256_cmpeq_epi64(aaaa2 & mask22, mask22);
            __m256i cmask52 = _mm256_cmpeq_epi64(aaaa2 & mask31, mask31);
            __m256i cmask62 = _mm256_cmpeq_epi64(aaaa2 & mask32, mask32);
            __m256i cmask72 = _mm256_cmpeq_epi64(aaaa2 & mask41, mask41);
            __m256i cmask82 = _mm256_cmpeq_epi64(aaaa2 & mask42, mask42);

            for (int r = 0; r <= c; r++) {

                int mp = MAYO_POS;

                __m256i in0swap = _mm256_permute4x64_epi64(P1[mp+0], 0b01001110);
                __m256i in1swap = _mm256_permute4x64_epi64(P1[mp+1], 0b01001110);

                acc[2*(r * K_MAX + k) + 0] ^= (P1[mp+0] & cmask1) ^ (in0swap  & cmask3) 
                                            ^ (P1[mp+1] & cmask5) ^ (in1swap  & cmask7);
                acc[2*(r * K_MAX + k) + 1] ^= (P1[mp+1] & cmask2) ^ (in1swap  & cmask4)
                                            ^ (P1[mp+0] & cmask6) ^ (in0swap  & cmask8);

                acc[2*(r * K_MAX + k + 1) + 0] ^= (P1[mp+0] & cmask12) ^ (in0swap  & cmask32) 
                                                ^ (P1[mp+1] & cmask52) ^ (in1swap  & cmask72);
                acc[2*(r * K_MAX + k + 1) + 1] ^= (P1[mp+1] & cmask22) ^ (in1swap  & cmask42)
                                                ^ (P1[mp+0] & cmask62) ^ (in0swap  & cmask82);
            }
        }
    }

    #undef MAYO_POS
}


#endif