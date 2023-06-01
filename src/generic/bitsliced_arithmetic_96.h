// SPDX-License-Identifier: Apache-2.0

#ifndef BITSLICED_ARITHMETIC_96_H
#define BITSLICED_ARITHMETIC_96_H

#include <stdint.h>
#include <mayo.h>

// This implements arithmetic for bitsliced vectors of 96 field elements in Z_2[x]/(x^4+x+1)

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
inline void bitsliced_96_vec_sum(const uint32_t *in1, const uint32_t *in2, uint32_t *out) {
    out[0]  = in1[0] ^ in2[0];
    out[1]  = in1[1] ^ in2[1];
    out[2]  = in1[2] ^ in2[2];
    out[3]  = in1[3] ^ in2[3];
    out[4]  = in1[4] ^ in2[4];
    out[5]  = in1[5] ^ in2[5];
    out[6]  = in1[6] ^ in2[6];
    out[7]  = in1[7] ^ in2[7];
    out[8]  = in1[8] ^ in2[8];
    out[9]  = in1[9] ^ in2[9];
    out[10] = in1[10] ^ in2[10];
    out[11] = in1[11] ^ in2[11];
}

static
inline void bitsliced_96_vec_add(const uint32_t *in, uint32_t *acc) {
    acc[0] ^= in[0];
    acc[1] ^= in[1];
    acc[2] ^= in[2];
    acc[3] ^= in[3];
    acc[4] ^= in[4];
    acc[5] ^= in[5];
    acc[6] ^= in[6];
    acc[7] ^= in[7];
    acc[8] ^= in[8];
    acc[9] ^= in[9];
    acc[10] ^= in[10];
    acc[11] ^= in[11];
}

static 
inline void bitsliced_96_vec_mul_add_x(const uint32_t *in, uint32_t *acc) {
    acc[0]  ^= in[9];
    acc[1]  ^= in[10];
    acc[2]  ^= in[11];
    acc[3]  ^= in[0] ^ in[9];

    acc[4]  ^= in[1] ^ in[10];
    acc[5]  ^= in[2] ^ in[11];
    acc[6]  ^= in[3];
    acc[7]  ^= in[4];
    
    acc[8]  ^= in[5];
    acc[9]  ^= in[6];
    acc[10] ^= in[7];
    acc[11] ^= in[8];
}

static 
inline void bitsliced_96_vec_mul_add_x_inv(const uint32_t *in, uint32_t *acc) {
    acc[0]  ^= in[0] ^ in[3];
    acc[1]  ^= in[1] ^ in[4];
    acc[2]  ^= in[2] ^ in[5];
    acc[3]  ^= in[6];
    acc[4]  ^= in[7];
    acc[5]  ^= in[8];
    acc[6]  ^= in[9];
    acc[7]  ^= in[10];
    acc[8]  ^= in[11];
    acc[9]  ^= in[0];
    acc[10] ^= in[1];
    acc[11] ^= in[2];
}

static 
inline void bitsliced_96_vec_mul_add(const uint32_t *in, unsigned char a, uint32_t *acc) {
    // terms of a
    uint32_t aa = a;
    uint32_t a0 = -(aa & 1);
    uint32_t a1 = -((aa >> 1) & 1);
    uint32_t a2 = -((aa >> 2) & 1);
    uint32_t a3 = -((aa >> 3) & 1);

    uint64_t x[3], y[3], z[3];
    // deg 0 term of a;
    acc[0]  ^= a0 & in[0];
    acc[1]  ^= a0 & in[1];
    acc[2]  ^= a0 & in[2];
    acc[3]  ^= a0 & in[3];
    acc[4]  ^= a0 & in[4];
    acc[5]  ^= a0 & in[5];
    acc[6]  ^= a0 & in[6];
    acc[7]  ^= a0 & in[7];
    acc[8]  ^= a0 & in[8];
    acc[9]  ^= a0 & in[9];
    acc[10] ^= a0 & in[10];
    acc[11] ^= a0 & in[11];

    // deg 1 term of a;
    x[0] = in[0] ^ in[9];
    x[1] = in[1] ^ in[10];
    x[2] = in[2] ^ in[11];

    acc[0]  ^= a1 & in[9];
    acc[1]  ^= a1 & in[10];
    acc[2]  ^= a1 & in[11];
    acc[3]  ^= a1 & x[0];
    acc[4]  ^= a1 & x[1];
    acc[5]  ^= a1 & x[2];
    acc[6]  ^= a1 & in[3];
    acc[7]  ^= a1 & in[4];
    acc[8]  ^= a1 & in[5];
    acc[9]  ^= a1 & in[6];
    acc[10] ^= a1 & in[7];
    acc[11] ^= a1 & in[8];

    // deg 2 term of a
    y[0] = in[9] ^ in[6];
    y[1] = in[10] ^ in[7];
    y[2] = in[11] ^ in[8];

    acc[0]  ^= a2 & in[6];
    acc[1]  ^= a2 & in[7];
    acc[2]  ^= a2 & in[8];
    acc[3]  ^= a2 & y[0];
    acc[4]  ^= a2 & y[1];
    acc[5]  ^= a2 & y[2];
    acc[6]  ^= a2 & x[0];
    acc[7]  ^= a2 & x[1];
    acc[8]  ^= a2 & x[2];
    acc[9]  ^= a2 & in[3];
    acc[10] ^= a2 & in[4];
    acc[11] ^= a2 & in[5];

    // deg 3 term of a
    z[0] = in[6] ^ in[3];
    z[1] = in[7] ^ in[4];
    z[2] = in[8] ^ in[5];
    acc[0]  ^= a3 & in[3];
    acc[1]  ^= a3 & in[4];
    acc[2]  ^= a3 & in[5];
    acc[3]  ^= a3 & z[0];
    acc[4]  ^= a3 & z[1];
    acc[5]  ^= a3 & z[2];
    acc[6]  ^= a3 & y[0];
    acc[7]  ^= a3 & y[1];
    acc[8]  ^= a3 & y[2];
    acc[9]  ^= a3 & x[0];
    acc[10] ^= a3 & x[1];
    acc[11] ^= a3 & x[2];
}

static 
inline void bitsliced_96_vec_mul_add_a(const uint32_t *in, uint32_t a0, uint32_t a1, uint32_t a2, uint32_t a3, uint32_t *acc) {
    // terms of a
    uint32_t x[3], y[3], z[3];
    // deg 0 term of a;
    acc[0]  ^= a0 & in[0];
    acc[1]  ^= a0 & in[1];
    acc[2]  ^= a0 & in[2];
    acc[3]  ^= a0 & in[3];
    acc[4]  ^= a0 & in[4];
    acc[5]  ^= a0 & in[5];
    acc[6]  ^= a0 & in[6];
    acc[7]  ^= a0 & in[7];
    acc[8]  ^= a0 & in[8];
    acc[9]  ^= a0 & in[9];
    acc[10] ^= a0 & in[10];
    acc[11] ^= a0 & in[11];

    // deg 1 term of a;
    x[0] = in[0] ^ in[9];
    x[1] = in[1] ^ in[10];
    x[2] = in[2] ^ in[11];

    acc[0]  ^= a1 & in[9];
    acc[1]  ^= a1 & in[10];
    acc[2]  ^= a1 & in[11];
    acc[3]  ^= a1 & x[0];
    acc[4]  ^= a1 & x[1];
    acc[5]  ^= a1 & x[2];
    acc[6]  ^= a1 & in[3];
    acc[7]  ^= a1 & in[4];
    acc[8]  ^= a1 & in[5];
    acc[9]  ^= a1 & in[6];
    acc[10] ^= a1 & in[7];
    acc[11] ^= a1 & in[8];

    // deg 2 term of a
    y[0] = in[9] ^ in[6];
    y[1] = in[10] ^ in[7];
    y[2] = in[11] ^ in[8];

    acc[0]  ^= a2 & in[6];
    acc[1]  ^= a2 & in[7];
    acc[2]  ^= a2 & in[8];
    acc[3]  ^= a2 & y[0];
    acc[4]  ^= a2 & y[1];
    acc[5]  ^= a2 & y[2];
    acc[6]  ^= a2 & x[0];
    acc[7]  ^= a2 & x[1];
    acc[8]  ^= a2 & x[2];
    acc[9]  ^= a2 & in[3];
    acc[10] ^= a2 & in[4];
    acc[11] ^= a2 & in[5];

    // deg 3 term of a
    z[0] = in[6] ^ in[3];
    z[1] = in[7] ^ in[4];
    z[2] = in[8] ^ in[5];
    acc[0]  ^= a3 & in[3];
    acc[1]  ^= a3 & in[4];
    acc[2]  ^= a3 & in[5];
    acc[3]  ^= a3 & z[0];
    acc[4]  ^= a3 & z[1];
    acc[5]  ^= a3 & z[2];
    acc[6]  ^= a3 & y[0];
    acc[7]  ^= a3 & y[1];
    acc[8]  ^= a3 & y[2];
    acc[9]  ^= a3 & x[0];
    acc[10] ^= a3 & x[1];
    acc[11] ^= a3 & x[2];
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