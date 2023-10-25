// SPDX-License-Identifier: Apache-2.0

#ifndef BITSLICED_ARITHMETIC_64_H
#define BITSLICED_ARITHMETIC_64_H

#include <stdint.h>
#include <mayo.h>
#include <immintrin.h>

// This implements arithmetic for bitsliced vectors of 64 field elements in Z_2[x]/(x^4+x+1)

static
inline void bitsliced_64_vec_copy(const uint64_t *in, uint64_t *out) {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
    out[3] = in[3];
}

static
inline void bitsliced_64_vec_sum(const uint64_t *_in1, const uint64_t *_in2, uint64_t *_out) {
    const __m256i* in1 = (const __m256i*) _in1;
    const __m256i* in2 = (const __m256i*) _in2;
    __m256i* out = (__m256i*) _out;
    out[0] = in1[0] ^ in2[0];
}

static
inline void bitsliced_64_vec_add(const uint64_t *_in, uint64_t *_acc) {
    const __m256i* in = (const __m256i*) _in;
    __m256i* acc = (__m256i*)_acc;
    acc[0] ^= in[0];
}

static 
inline void bitsliced_64_vec_mul_add_x(const uint64_t* _in, uint64_t* _acc) {
    const __m256i* in = (const __m256i*) _in;
    __m256i* acc = (__m256i*)_acc;

    const __m256i andmask = _mm256_set_epi64x(0, 0, 0xffffffffffffffff, 0);
    acc[0] ^= _mm256_permute4x64_epi64(in[0], 0b10011111) ^ (_mm256_shuffle_epi32(in[0], 0b01001110) & andmask);
}

static 
inline void bitsliced_64_vec_mul_add_x_inv(const uint64_t *_in, uint64_t *_acc) {
    const __m256i* in = (const __m256i*) _in;
    __m256i* acc = (__m256i*)_acc;

    const __m256i andmask = _mm256_set_epi64x(0, 0, 0, 0xffffffffffffffff);
    acc[0] ^= _mm256_permute4x64_epi64(in[0], 0b00111001) ^ (in[0] & andmask);
}

static
inline void bitsliced_64_vec_mul_add(const uint64_t *_in, unsigned char a, uint64_t *_acc) {
    const __m256i* in = (const __m256i*) _in;
    __m256i* acc = (__m256i*) _acc;

    const __m256i lut_a = _mm256_setr_epi8(
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1,
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1
    );


    __m256i inshuf2 = _mm256_permute4x64_epi64(in[0], 0b01001110);
    __m256i inshuf1 = _mm256_shuffle_epi32(in[0], 0b01001110);
    __m256i inshuf3 = _mm256_shuffle_epi32(inshuf2, 0b01001110); 

    const __m256i mask1 = _mm256_set_epi64x(16,  16, 16, 1  );
    const __m256i mask2 = _mm256_set_epi64x(32,  8,  32, 128);
    const __m256i mask3 = _mm256_set_epi64x(64,  64, 4,  64 );
    const __m256i mask4 = _mm256_set_epi64x(128, 32, 8,  32 );

    __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(a));

    acc[0] ^= in[0]   & _mm256_cmpeq_epi64(aaaa & mask1, mask1);
    acc[0] ^= inshuf1 & _mm256_cmpeq_epi64(aaaa & mask2, mask2);
    acc[0] ^= inshuf2 & _mm256_cmpeq_epi64(aaaa & mask3, mask3);
    acc[0] ^= inshuf3 & _mm256_cmpeq_epi64(aaaa & mask4, mask4);
}

static
inline void mayo_12_P1_times_O(const uint32_t *_P1, const unsigned char *O, uint32_t *_acc) {

    const __m256i* P1 = (const __m256i*) _P1;
    __m256i* acc = (__m256i*) _acc;

    const __m256i lut_a = _mm256_setr_epi8(
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1,
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1
    );
    const __m256i mask1 = _mm256_set_epi64x(16,  16, 16, 1  );
    const __m256i mask2 = _mm256_set_epi64x(32,  8,  32, 128);
    const __m256i mask3 = _mm256_set_epi64x(64,  64, 4,  64 );
    const __m256i mask4 = _mm256_set_epi64x(128, 32, 8,  32 );

    #define MAYO_POS ( r*V_MAX + c - (r)*(r+1)/2 )

    for (int c = 0; c < V_MAX; c++) {
        int k;
        for (k = 0; k < (O_MAX/4)*4; k += 4) {

            __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k]));
            __m256i cmask1 = _mm256_cmpeq_epi64(aaaa & mask1, mask1);
            __m256i cmask2 = _mm256_cmpeq_epi64(aaaa & mask2, mask2);
            __m256i cmask3 = _mm256_cmpeq_epi64(aaaa & mask3, mask3);
            __m256i cmask4 = _mm256_cmpeq_epi64(aaaa & mask4, mask4);

            __m256i aaaa2 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 1]));
            __m256i cmask12 = _mm256_cmpeq_epi64(aaaa2 & mask1, mask1);
            __m256i cmask22 = _mm256_cmpeq_epi64(aaaa2 & mask2, mask2);
            __m256i cmask32 = _mm256_cmpeq_epi64(aaaa2 & mask3, mask3);
            __m256i cmask42 = _mm256_cmpeq_epi64(aaaa2 & mask4, mask4);

            __m256i aaaa3 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 2]));
            __m256i cmask13 = _mm256_cmpeq_epi64(aaaa3 & mask1, mask1);
            __m256i cmask23 = _mm256_cmpeq_epi64(aaaa3 & mask2, mask2);
            __m256i cmask33 = _mm256_cmpeq_epi64(aaaa3 & mask3, mask3);
            __m256i cmask43 = _mm256_cmpeq_epi64(aaaa3 & mask4, mask4);

            __m256i aaaa4 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 3]));
            __m256i cmask14 = _mm256_cmpeq_epi64(aaaa4 & mask1, mask1);
            __m256i cmask24 = _mm256_cmpeq_epi64(aaaa4 & mask2, mask2);
            __m256i cmask34 = _mm256_cmpeq_epi64(aaaa4 & mask3, mask3);
            __m256i cmask44 = _mm256_cmpeq_epi64(aaaa4 & mask4, mask4);

            for (int r = 0; r <= c; r++) {

                int mp = MAYO_POS;

                __m256i inshuf2 = _mm256_permute4x64_epi64(P1[mp], 0b01001110);
                __m256i inshuf1 = _mm256_shuffle_epi32(P1[mp], 0b01001110);
                __m256i inshuf3 = _mm256_shuffle_epi32(inshuf2, 0b01001110); 
                
                acc[r * O_MAX + k] ^= P1[mp]  & cmask1;
                acc[r * O_MAX + k] ^= inshuf1 & cmask2;
                acc[r * O_MAX + k] ^= inshuf2 & cmask3;
                acc[r * O_MAX + k] ^= inshuf3 & cmask4;

                acc[r * O_MAX + k + 1] ^= P1[mp]  & cmask12;
                acc[r * O_MAX + k + 1] ^= inshuf1 & cmask22;
                acc[r * O_MAX + k + 1] ^= inshuf2 & cmask32;
                acc[r * O_MAX + k + 1] ^= inshuf3 & cmask42;

                acc[r * O_MAX + k + 2] ^= P1[mp] & cmask13;
                acc[r * O_MAX + k + 2] ^= inshuf1            & cmask23;
                acc[r * O_MAX + k + 2] ^= inshuf2            & cmask33;
                acc[r * O_MAX + k + 2] ^= inshuf3            & cmask43;

                acc[r * O_MAX + k + 3] ^= P1[mp] & cmask14;
                acc[r * O_MAX + k + 3] ^= inshuf1            & cmask24;
                acc[r * O_MAX + k + 3] ^= inshuf2            & cmask34;
                acc[r * O_MAX + k + 3] ^= inshuf3            & cmask44;
            }
        }
        for (; k < (O_MAX/2)*2; k += 2) { // if o is 2 mod 4

            __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k]));
            __m256i cmask1 = _mm256_cmpeq_epi64(aaaa & mask1, mask1);
            __m256i cmask2 = _mm256_cmpeq_epi64(aaaa & mask2, mask2);
            __m256i cmask3 = _mm256_cmpeq_epi64(aaaa & mask3, mask3);
            __m256i cmask4 = _mm256_cmpeq_epi64(aaaa & mask4, mask4);

            __m256i aaaa2 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 1]));
            __m256i cmask12 = _mm256_cmpeq_epi64(aaaa2 & mask1, mask1);
            __m256i cmask22 = _mm256_cmpeq_epi64(aaaa2 & mask2, mask2);
            __m256i cmask32 = _mm256_cmpeq_epi64(aaaa2 & mask3, mask3);
            __m256i cmask42 = _mm256_cmpeq_epi64(aaaa2 & mask4, mask4);

            for (int r = 0; r <= c; r++) {

                int mp = MAYO_POS;

                __m256i inshuf2 = _mm256_permute4x64_epi64(P1[mp], 0b01001110);
                __m256i inshuf1 = _mm256_shuffle_epi32(P1[mp], 0b01001110);
                __m256i inshuf3 = _mm256_shuffle_epi32(inshuf2, 0b01001110); 
                
                acc[r * O_MAX + k] ^= P1[mp]  & cmask1;
                acc[r * O_MAX + k] ^= inshuf1 & cmask2;
                acc[r * O_MAX + k] ^= inshuf2 & cmask3;
                acc[r * O_MAX + k] ^= inshuf3 & cmask4;

                acc[r * O_MAX + k + 1] ^= P1[mp]  & cmask12;
                acc[r * O_MAX + k + 1] ^= inshuf1 & cmask22;
                acc[r * O_MAX + k + 1] ^= inshuf2 & cmask32;
                acc[r * O_MAX + k + 1] ^= inshuf3 & cmask42;
            }
        }
    }

    #undef MAYO_POS
}

static
inline void mayo_12_P1P1t_times_O(const uint32_t *_P1P1t, const unsigned char *O, uint32_t *_acc) {
    const __m256i* P1P1t = (const __m256i*) _P1P1t;
    __m256i* acc = (__m256i*) _acc;

    const __m256i lut_a = _mm256_setr_epi8(
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1,
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1
    );
    const __m256i mask1 = _mm256_set_epi64x(16,  16, 16, 1  );
    const __m256i mask2 = _mm256_set_epi64x(32,  8,  32, 128);
    const __m256i mask3 = _mm256_set_epi64x(64,  64, 4,  64 );
    const __m256i mask4 = _mm256_set_epi64x(128, 32, 8,  32 );
    
    for (int c = 0; c < V_MAX; c++) {
        int k;
        for (k = 0; k < (O_MAX/4)*4; k += 4) {
            __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k]));
            __m256i cmask1 = _mm256_cmpeq_epi64(aaaa & mask1, mask1);
            __m256i cmask2 = _mm256_cmpeq_epi64(aaaa & mask2, mask2);
            __m256i cmask3 = _mm256_cmpeq_epi64(aaaa & mask3, mask3);
            __m256i cmask4 = _mm256_cmpeq_epi64(aaaa & mask4, mask4);

            __m256i aaaa2 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 1]));
            __m256i cmask12 = _mm256_cmpeq_epi64(aaaa2 & mask1, mask1);
            __m256i cmask22 = _mm256_cmpeq_epi64(aaaa2 & mask2, mask2);
            __m256i cmask32 = _mm256_cmpeq_epi64(aaaa2 & mask3, mask3);
            __m256i cmask42 = _mm256_cmpeq_epi64(aaaa2 & mask4, mask4);

            __m256i aaaa3 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 2]));
            __m256i cmask13 = _mm256_cmpeq_epi64(aaaa3 & mask1, mask1);
            __m256i cmask23 = _mm256_cmpeq_epi64(aaaa3 & mask2, mask2);
            __m256i cmask33 = _mm256_cmpeq_epi64(aaaa3 & mask3, mask3);
            __m256i cmask43 = _mm256_cmpeq_epi64(aaaa3 & mask4, mask4);

            __m256i aaaa4 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 3]));
            __m256i cmask14 = _mm256_cmpeq_epi64(aaaa4 & mask1, mask1);
            __m256i cmask24 = _mm256_cmpeq_epi64(aaaa4 & mask2, mask2);
            __m256i cmask34 = _mm256_cmpeq_epi64(aaaa4 & mask3, mask3);
            __m256i cmask44 = _mm256_cmpeq_epi64(aaaa4 & mask4, mask4);

            for (int r = 0; r < V_MAX; r++) {

                __m256i inshuf2 = _mm256_permute4x64_epi64(P1P1t[r*V_MAX + c], 0b01001110);
                __m256i inshuf1 = _mm256_shuffle_epi32(P1P1t[r*V_MAX + c], 0b01001110);
                __m256i inshuf3 = _mm256_shuffle_epi32(inshuf2, 0b01001110); 

                __m256i acc0 = _mm256_loadu_si256(&acc[r * O_MAX + k]);
                __m256i acc1 = _mm256_loadu_si256(&acc[r * O_MAX + k + 1]);
                __m256i acc2 = _mm256_loadu_si256(&acc[r * O_MAX + k + 2]);
                __m256i acc3 = _mm256_loadu_si256(&acc[r * O_MAX + k + 3]);

                _mm256_storeu_si256(&acc[r * O_MAX + k],
                    acc0 ^ 
                    (P1P1t[r*V_MAX + c] & cmask1) ^
                    (inshuf1            & cmask2) ^
                    (inshuf2            & cmask3) ^
                    (inshuf3            & cmask4)
                );
                
                _mm256_storeu_si256(&acc[r * O_MAX + k + 1],
                    acc1 ^ 
                    (P1P1t[r*V_MAX + c] & cmask12) ^
                    (inshuf1            & cmask22) ^
                    (inshuf2            & cmask32) ^
                    (inshuf3            & cmask42)
                );

                _mm256_storeu_si256(&acc[r * O_MAX + k + 2],
                    acc2 ^ 
                    (P1P1t[r*V_MAX + c] & cmask13) ^
                    (inshuf1            & cmask23) ^
                    (inshuf2            & cmask33) ^
                    (inshuf3            & cmask43)
                );

                _mm256_storeu_si256(&acc[r * O_MAX + k + 3],
                    acc3 ^ 
                    (P1P1t[r*V_MAX + c] & cmask14) ^
                    (inshuf1            & cmask24) ^
                    (inshuf2            & cmask34) ^
                    (inshuf3            & cmask44)
                );
            }
        }
        for (; k < (O_MAX/2)*2; k += 2) {
            __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k]));
            __m256i cmask1 = _mm256_cmpeq_epi64(aaaa & mask1, mask1);
            __m256i cmask2 = _mm256_cmpeq_epi64(aaaa & mask2, mask2);
            __m256i cmask3 = _mm256_cmpeq_epi64(aaaa & mask3, mask3);
            __m256i cmask4 = _mm256_cmpeq_epi64(aaaa & mask4, mask4);

            __m256i aaaa2 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(O[c * O_MAX + k + 1]));
            __m256i cmask12 = _mm256_cmpeq_epi64(aaaa2 & mask1, mask1);
            __m256i cmask22 = _mm256_cmpeq_epi64(aaaa2 & mask2, mask2);
            __m256i cmask32 = _mm256_cmpeq_epi64(aaaa2 & mask3, mask3);
            __m256i cmask42 = _mm256_cmpeq_epi64(aaaa2 & mask4, mask4);

            for (int r = 0; r < V_MAX; r++) {

                __m256i inshuf2 = _mm256_permute4x64_epi64(P1P1t[r*V_MAX + c], 0b01001110);
                __m256i inshuf1 = _mm256_shuffle_epi32(P1P1t[r*V_MAX + c], 0b01001110);
                __m256i inshuf3 = _mm256_shuffle_epi32(inshuf2, 0b01001110); 

                __m256i acc0 = _mm256_loadu_si256(&acc[r * O_MAX + k]);
                __m256i acc1 = _mm256_loadu_si256(&acc[r * O_MAX + k + 1]);

                _mm256_storeu_si256(&acc[r * O_MAX + k],
                    acc0 ^ 
                    (P1P1t[r*V_MAX + c] & cmask1) ^
                    (inshuf1            & cmask2) ^
                    (inshuf2            & cmask3) ^
                    (inshuf3            & cmask4)
                );

                _mm256_storeu_si256(&acc[r * O_MAX + k + 1],
                    acc1 ^ 
                    (P1P1t[r*V_MAX + c] & cmask12) ^
                    (inshuf1            & cmask22) ^
                    (inshuf2            & cmask32) ^
                    (inshuf3            & cmask42)
                );
            }
        }
    }
}

// multiplies m bitsliced (possibly upper triangular) matrices with the transpose of a single matrix and adds result to acc
static
inline void mayo_12_P1_times_Vt(const uint32_t *_P1, const unsigned char *V, uint32_t *_acc) {
    const __m256i* P1 = (const __m256i*) _P1;
    __m256i* acc = (__m256i*) _acc;

    const __m256i lut_a = _mm256_setr_epi8(
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1,
        0x0, 0x13, 0x26, 0x35, 0x4c, 0x5f, 0x6a, 0x79, 0x98, 0x8b, 0xbe, 0xad, 0xd4, 0xc7, 0xf2, 0xe1
    );
    const __m256i mask1 = _mm256_set_epi64x(16,  16, 16, 1  );
    const __m256i mask2 = _mm256_set_epi64x(32,  8,  32, 128);
    const __m256i mask3 = _mm256_set_epi64x(64,  64, 4,  64 );
    const __m256i mask4 = _mm256_set_epi64x(128, 32, 8,  32 );

    #define MAYO_POS ( r*V_MAX + c - (r)*(r+1)/2 )

    for (int c = 0; c < V_MAX; c++) {
        
        __m256i aaaa = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[c]));
        __m256i cmask1 = _mm256_cmpeq_epi64(aaaa & mask1, mask1);
        __m256i cmask2 = _mm256_cmpeq_epi64(aaaa & mask2, mask2);
        __m256i cmask3 = _mm256_cmpeq_epi64(aaaa & mask3, mask3);
        __m256i cmask4 = _mm256_cmpeq_epi64(aaaa & mask4, mask4);

        __m256i aaaa2 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[V_MAX + c ]));
        __m256i cmask12 = _mm256_cmpeq_epi64(aaaa2 & mask1, mask1);
        __m256i cmask22 = _mm256_cmpeq_epi64(aaaa2 & mask2, mask2);
        __m256i cmask32 = _mm256_cmpeq_epi64(aaaa2 & mask3, mask3);
        __m256i cmask42 = _mm256_cmpeq_epi64(aaaa2 & mask4, mask4);

        __m256i aaaa3 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[2 * V_MAX + c]));
        __m256i cmask13 = _mm256_cmpeq_epi64(aaaa3 & mask1, mask1);
        __m256i cmask23 = _mm256_cmpeq_epi64(aaaa3 & mask2, mask2);
        __m256i cmask33 = _mm256_cmpeq_epi64(aaaa3 & mask3, mask3);
        __m256i cmask43 = _mm256_cmpeq_epi64(aaaa3 & mask4, mask4);

        __m256i aaaa4 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[3 * V_MAX + c]));
        __m256i cmask14 = _mm256_cmpeq_epi64(aaaa4 & mask1, mask1);
        __m256i cmask24 = _mm256_cmpeq_epi64(aaaa4 & mask2, mask2);
        __m256i cmask34 = _mm256_cmpeq_epi64(aaaa4 & mask3, mask3);
        __m256i cmask44 = _mm256_cmpeq_epi64(aaaa4 & mask4, mask4);

        for (int r = 0; r <= c; r++) {

            int mp = MAYO_POS;

            __m256i inshuf2 = _mm256_permute4x64_epi64(P1[mp], 0b01001110);
            __m256i inshuf1 = _mm256_shuffle_epi32(P1[mp], 0b01001110);
            __m256i inshuf3 = _mm256_shuffle_epi32(inshuf2, 0b01001110); 
            
            acc[r * K_MAX] ^= P1[mp]  & cmask1;
            acc[r * K_MAX] ^= inshuf1 & cmask2;
            acc[r * K_MAX] ^= inshuf2 & cmask3;
            acc[r * K_MAX] ^= inshuf3 & cmask4;

            acc[r * K_MAX + 1] ^= P1[mp]  & cmask12;
            acc[r * K_MAX + 1] ^= inshuf1 & cmask22;
            acc[r * K_MAX + 1] ^= inshuf2 & cmask32;
            acc[r * K_MAX + 1] ^= inshuf3 & cmask42;

            acc[r * K_MAX + 2] ^= P1[mp] & cmask13;
            acc[r * K_MAX + 2] ^= inshuf1            & cmask23;
            acc[r * K_MAX + 2] ^= inshuf2            & cmask33;
            acc[r * K_MAX + 2] ^= inshuf3            & cmask43;

            acc[r * K_MAX + 3] ^= P1[mp] & cmask14;
            acc[r * K_MAX + 3] ^= inshuf1            & cmask24;
            acc[r * K_MAX + 3] ^= inshuf2            & cmask34;
            acc[r * K_MAX + 3] ^= inshuf3            & cmask44;
        }
#if (K_MAX == 9)
        __m256i aaaa5 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[4 * V_MAX + c]));
        __m256i cmask15 = _mm256_cmpeq_epi64(aaaa5 & mask1, mask1);
        __m256i cmask25 = _mm256_cmpeq_epi64(aaaa5 & mask2, mask2);
        __m256i cmask35 = _mm256_cmpeq_epi64(aaaa5 & mask3, mask3);
        __m256i cmask45 = _mm256_cmpeq_epi64(aaaa5 & mask4, mask4);

        __m256i aaaa6 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[5 * V_MAX + c]));
        __m256i cmask16 = _mm256_cmpeq_epi64(aaaa6 & mask1, mask1);
        __m256i cmask26 = _mm256_cmpeq_epi64(aaaa6 & mask2, mask2);
        __m256i cmask36 = _mm256_cmpeq_epi64(aaaa6 & mask3, mask3);
        __m256i cmask46 = _mm256_cmpeq_epi64(aaaa6 & mask4, mask4);

        __m256i aaaa7 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[6 * V_MAX + c]));
        __m256i cmask17 = _mm256_cmpeq_epi64(aaaa7 & mask1, mask1);
        __m256i cmask27 = _mm256_cmpeq_epi64(aaaa7 & mask2, mask2);
        __m256i cmask37 = _mm256_cmpeq_epi64(aaaa7 & mask3, mask3);
        __m256i cmask47 = _mm256_cmpeq_epi64(aaaa7 & mask4, mask4);

        __m256i aaaa8 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[7 * V_MAX + c]));
        __m256i cmask18 = _mm256_cmpeq_epi64(aaaa8 & mask1, mask1);
        __m256i cmask28 = _mm256_cmpeq_epi64(aaaa8 & mask2, mask2);
        __m256i cmask38 = _mm256_cmpeq_epi64(aaaa8 & mask3, mask3);
        __m256i cmask48 = _mm256_cmpeq_epi64(aaaa8 & mask4, mask4);

        __m256i aaaa9 = _mm256_shuffle_epi8(lut_a, _mm256_set1_epi8(V[8 * V_MAX + c]));
        __m256i cmask19 = _mm256_cmpeq_epi64(aaaa9 & mask1, mask1);
        __m256i cmask29 = _mm256_cmpeq_epi64(aaaa9 & mask2, mask2);
        __m256i cmask39 = _mm256_cmpeq_epi64(aaaa9 & mask3, mask3);
        __m256i cmask49 = _mm256_cmpeq_epi64(aaaa9 & mask4, mask4);

        for (int r = 0; r <= c; r++) {

            int mp = MAYO_POS;

            __m256i inshuf2 = _mm256_permute4x64_epi64(P1[mp], 0b01001110);
            __m256i inshuf1 = _mm256_shuffle_epi32(P1[mp], 0b01001110);
            __m256i inshuf3 = _mm256_shuffle_epi32(inshuf2, 0b01001110); 
            
            acc[r * K_MAX + 4] ^= P1[mp]  & cmask15;
            acc[r * K_MAX + 4] ^= inshuf1 & cmask25;
            acc[r * K_MAX + 4] ^= inshuf2 & cmask35;
            acc[r * K_MAX + 4] ^= inshuf3 & cmask45;

            acc[r * K_MAX + 5] ^= P1[mp]  & cmask16;
            acc[r * K_MAX + 5] ^= inshuf1 & cmask26;
            acc[r * K_MAX + 5] ^= inshuf2 & cmask36;
            acc[r * K_MAX + 5] ^= inshuf3 & cmask46;

            acc[r * K_MAX + 6] ^= P1[mp]  & cmask17;
            acc[r * K_MAX + 6] ^= inshuf1 & cmask27;
            acc[r * K_MAX + 6] ^= inshuf2 & cmask37;
            acc[r * K_MAX + 6] ^= inshuf3 & cmask47;

            acc[r * K_MAX + 7] ^= P1[mp]  & cmask18;
            acc[r * K_MAX + 7] ^= inshuf1 & cmask28;
            acc[r * K_MAX + 7] ^= inshuf2 & cmask38;
            acc[r * K_MAX + 7] ^= inshuf3 & cmask48;

            acc[r * K_MAX + 8] ^= P1[mp]  & cmask19;
            acc[r * K_MAX + 8] ^= inshuf1 & cmask29;
            acc[r * K_MAX + 8] ^= inshuf2 & cmask39;
            acc[r * K_MAX + 8] ^= inshuf3 & cmask49;
        }
#endif
    }
    #undef MAYO_POS
}

static 
inline void bitsliced_64_multiply_bins(uint32_t *bins_32, uint32_t *out_32) {

    uint64_t *bins = (uint64_t *) bins_32;
    uint64_t *out = (uint64_t *) out_32;

    bitsliced_64_vec_mul_add_x_inv(bins +  5 * 4, bins +  10 * 4);
    bitsliced_64_vec_mul_add_x(bins + 11 * 4, bins + 12 * 4);
    bitsliced_64_vec_mul_add_x_inv(bins +  10 * 4, bins +  7 * 4);
    bitsliced_64_vec_mul_add_x(bins + 12 * 4, bins +  6 * 4);
    bitsliced_64_vec_mul_add_x_inv(bins +  7 * 4, bins +  14 * 4);
    bitsliced_64_vec_mul_add_x(bins +  6 * 4, bins +  3 * 4);
    bitsliced_64_vec_mul_add_x_inv(bins +  14 * 4, bins +  15 * 4);
    bitsliced_64_vec_mul_add_x(bins +  3 * 4, bins +  8 * 4);
    bitsliced_64_vec_mul_add_x_inv(bins +  15 * 4, bins +  13 * 4);
    bitsliced_64_vec_mul_add_x(bins +  8 * 4, bins +  4 * 4);
    bitsliced_64_vec_mul_add_x_inv(bins +  13 * 4, bins +  9 * 4);
    bitsliced_64_vec_mul_add_x(bins +  4 * 4, bins +  2 * 4);
    bitsliced_64_vec_mul_add_x_inv(bins +   9 * 4, bins +  1 * 4);
    bitsliced_64_vec_mul_add_x(bins +  2 * 4, bins +  1 * 4);
    bitsliced_64_vec_copy(bins + 4, out);
}

#endif