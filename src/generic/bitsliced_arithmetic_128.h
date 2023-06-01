// SPDX-License-Identifier: Apache-2.0

#ifndef BITSLICED_ARITHMETIC_128_H
#define BITSLICED_ARITHMETIC_128_H

#include <stdint.h>
#include <mayo.h>

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
inline void bitsliced_128_vec_sum(const uint64_t *in1, const uint64_t *in2, uint64_t *out) {
    out[0] = in1[0] ^ in2[0];
    out[1] = in1[1] ^ in2[1];
    out[2] = in1[2] ^ in2[2];
    out[3] = in1[3] ^ in2[3];
    out[4] = in1[4] ^ in2[4];
    out[5] = in1[5] ^ in2[5];
    out[6] = in1[6] ^ in2[6];
    out[7] = in1[7] ^ in2[7];
}

static
inline void bitsliced_128_vec_add(const uint64_t *in, uint64_t *acc) {
    acc[0] ^= in[0];
    acc[1] ^= in[1];
    acc[2] ^= in[2];
    acc[3] ^= in[3];
    acc[4] ^= in[4];
    acc[5] ^= in[5];
    acc[6] ^= in[6];
    acc[7] ^= in[7];
}

static 
inline void bitsliced_128_vec_mul_add_x(const uint64_t *in, uint64_t *acc) {
    acc[0] ^= in[6];
    acc[1] ^= in[7];
    acc[2] ^= in[0] ^ in[6];
    acc[3] ^= in[1] ^ in[7];
    acc[4] ^= in[2];
    acc[5] ^= in[3];
    acc[6] ^= in[4];
    acc[7] ^= in[5];
}

static 
inline void bitsliced_128_vec_mul_add_x_inv(const uint64_t *in, uint64_t *acc) {
    acc[0] ^= in[2] ^ in[0];
    acc[1] ^= in[3] ^ in[1];
    acc[2] ^= in[4];
    acc[3] ^= in[5];
    acc[4] ^= in[6];
    acc[5] ^= in[7];
    acc[6] ^= in[0];
    acc[7] ^= in[1];
}

static 
inline void bitsliced_128_vec_mul_add(const uint64_t *in, unsigned char a, uint64_t *acc) {
    // terms of a
    uint64_t aa = a;
    uint64_t a0 = -(aa & 1);
    uint64_t a1 = -((aa >> 1) & 1);
    uint64_t a2 = -((aa >> 2) & 1);
    uint64_t a3 = -((aa >> 3) & 1);

    uint64_t x[2], y[2], z[2];
    // deg 0 term of a;
    acc[0] ^= a0 & in[0];
    acc[1] ^= a0 & in[1];
    acc[2] ^= a0 & in[2];
    acc[3] ^= a0 & in[3];
    acc[4] ^= a0 & in[4];
    acc[5] ^= a0 & in[5];
    acc[6] ^= a0 & in[6];
    acc[7] ^= a0 & in[7];

    // deg 1 term of a;
    x[0] = in[0] ^ in[6];
    x[1] = in[1] ^ in[7];
    acc[0] ^= a1 & in[6];
    acc[1] ^= a1 & in[7];
    acc[2] ^= a1 & x[0];
    acc[3] ^= a1 & x[1];
    acc[4] ^= a1 & in[2];
    acc[5] ^= a1 & in[3];
    acc[6] ^= a1 & in[4];
    acc[7] ^= a1 & in[5];

    // deg 2 term of a
    y[0] = in[6] ^ in[4];
    y[1] = in[7] ^ in[5];
    acc[0] ^= a2 & in[4];
    acc[1] ^= a2 & in[5];
    acc[2] ^= a2 & y[0];
    acc[3] ^= a2 & y[1];
    acc[4] ^= a2 & x[0];
    acc[5] ^= a2 & x[1];
    acc[6] ^= a2 & in[2];
    acc[7] ^= a2 & in[3];

    // deg 3 term of a
    z[0] = in[4] ^ in[2];
    z[1] = in[5] ^ in[3];
    acc[0] ^= a3 & in[2];
    acc[1] ^= a3 & in[3];
    acc[2] ^= a3 & z[0];
    acc[3] ^= a3 & z[1];
    acc[4] ^= a3 & y[0];
    acc[5] ^= a3 & y[1];
    acc[6] ^= a3 & x[0];
    acc[7] ^= a3 & x[1];
}

static 
inline void bitsliced_128_vec_mul_add_a(const uint64_t *in, uint64_t a0, uint64_t a1, uint64_t a2, uint64_t a3, uint64_t *acc) {

    uint64_t x[2], y[2], z[2];
    // deg 0 term of a;
    acc[0] ^= a0 & in[0];
    acc[1] ^= a0 & in[1];
    acc[2] ^= a0 & in[2];
    acc[3] ^= a0 & in[3];
    acc[4] ^= a0 & in[4];
    acc[5] ^= a0 & in[5];
    acc[6] ^= a0 & in[6];
    acc[7] ^= a0 & in[7];

    // deg 1 term of a;
    x[0] = in[0] ^ in[6];
    x[1] = in[1] ^ in[7];
    acc[0] ^= a1 & in[6];
    acc[1] ^= a1 & in[7];
    acc[2] ^= a1 & x[0];
    acc[3] ^= a1 & x[1];
    acc[4] ^= a1 & in[2];
    acc[5] ^= a1 & in[3];
    acc[6] ^= a1 & in[4];
    acc[7] ^= a1 & in[5];

    // deg 2 term of a
    y[0] = in[6] ^ in[4];
    y[1] = in[7] ^ in[5];
    acc[0] ^= a2 & in[4];
    acc[1] ^= a2 & in[5];
    acc[2] ^= a2 & y[0];
    acc[3] ^= a2 & y[1];
    acc[4] ^= a2 & x[0];
    acc[5] ^= a2 & x[1];
    acc[6] ^= a2 & in[2];
    acc[7] ^= a2 & in[3];

    // deg 3 term of a
    z[0] = in[4] ^ in[2];
    z[1] = in[5] ^ in[3];
    acc[0] ^= a3 & in[2];
    acc[1] ^= a3 & in[3];
    acc[2] ^= a3 & z[0];
    acc[3] ^= a3 & z[1];
    acc[4] ^= a3 & y[0];
    acc[5] ^= a3 & y[1];
    acc[6] ^= a3 & x[0];
    acc[7] ^= a3 & x[1];
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

#endif