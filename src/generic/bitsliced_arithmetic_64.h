// SPDX-License-Identifier: Apache-2.0

#ifndef BITSLICED_ARITHMETIC_64_H
#define BITSLICED_ARITHMETIC_64_H

#include <stdint.h>
#include <mayo.h>

// This implements arithmetic for bitsliced vectors of 64 field elements in Z_2[x]/(x^4+x+1)

static
inline void bitsliced_64_vec_copy(const uint64_t *in, uint64_t *out) {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
    out[3] = in[3];
}

static
inline void bitsliced_64_vec_sum(const uint64_t *in1, const uint64_t *in2, uint64_t *out) {
    out[0] = in1[0] ^ in2[0];
    out[1] = in1[1] ^ in2[1];
    out[2] = in1[2] ^ in2[2];
    out[3] = in1[3] ^ in2[3];
}

static
inline void bitsliced_64_vec_add(const uint64_t *in, uint64_t *acc) {
    acc[0] ^= in[0];
    acc[1] ^= in[1];
    acc[2] ^= in[2];
    acc[3] ^= in[3];
}

static 
inline void bitsliced_64_vec_mul_add_x(const uint64_t *in, uint64_t *acc) {
    acc[0] ^= in[3];
    acc[1] ^= in[0] ^ in[3];
    acc[2] ^= in[1];
    acc[3] ^= in[2];
}

static 
inline void bitsliced_64_vec_mul_add_x_inv(const uint64_t *in, uint64_t *acc) {
    acc[0] ^= in[0] ^ in[1];
    acc[1] ^= in[2];
    acc[2] ^= in[3];
    acc[3] ^= in[0];
}

static 
inline void bitsliced_64_vec_mul_add(const uint64_t *in, unsigned char a, uint64_t *acc) {
    // terms of a
    uint64_t aa = a;
    uint64_t a0 = -(aa & 1);
    uint64_t a1 = -((aa >> 1) & 1);
    uint64_t a2 = -((aa >> 2) & 1);
    uint64_t a3 = -((aa >> 3) & 1);

    uint64_t x, y, z;
    // deg 0 term of a;
    acc[0] ^= a0 & in[0];
    acc[1] ^= a0 & in[1];
    acc[2] ^= a0 & in[2];
    acc[3] ^= a0 & in[3];

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