
// SPDX-License-Identifier: Apache-2.0

#ifndef BITSLICED_ARITHMETIC_H
#define BITSLICED_ARITHMETIC_H

#include <bitsliced_arithmetic_64.h>
#include <bitsliced_arithmetic_96.h>
#include <bitsliced_arithmetic_128.h>
#include <stdint.h>
#include <mayo.h>
#include <stdint.h>

// This implements arithmetic for bitsliced vectors of m field elements in
// Z_2[x]/(x^4+x+1)

static
inline void bitsliced_m_vec_copy(int m_legs, const uint32_t *in,
                                 uint32_t *out) {
    for (int i = 0; i < m_legs * 4; i++) {
        out[i] = in[i];
    }
}

static
inline void bitsliced_m_vec_sum(int m_legs, const uint32_t *in1,
                                const uint32_t *in2, uint32_t *out) {
    for (int i = 0; i < m_legs * 4; i++) {
        out[i] = in1[i] ^ in2[i];
    }
}

static
inline void bitsliced_m_vec_add(int m_legs, const uint32_t *in,
                                uint32_t *acc) {
    for (int i = 0; i < m_legs * 4; i++) {
        acc[i] ^= in[i];
    }
}

void bitsliced_m_upper(int m_legs, const uint32_t *in, uint32_t *out, int size);
void bitsliced_m_vec_mul_add(int m_legs, const uint32_t *in, unsigned char a, uint32_t *acc);
void bitslice_m_vec(int m_legs, const unsigned char *in, uint32_t *out);

void unbitslice_m_vec(int m_legs, const uint32_t *in, unsigned char *out);

void mul_add_bitsliced_m_upper_triangular_mat_x_mat(
    int m_legs, const uint32_t *bs_mat, const unsigned char *mat, uint32_t *acc,
    int bs_mat_rows, int bs_mat_cols, int mat_cols, int triangular);

void mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(
    int m_legs, const uint32_t *bs_mat, const unsigned char *mat, uint32_t *acc,
    int bs_mat_rows, int bs_mat_cols, int mat_rows, int triangular);

void mul_add_mat_trans_x_bitsliced_m_mat(int m_legs, const unsigned char *mat,
        const uint32_t *bs_mat, uint32_t *acc,
        int mat_rows, int mat_cols,
        int bs_mat_cols);

void mul_add_mat_x_bitsliced_m_mat(int m_legs, const unsigned char *mat,
                                   const uint32_t *bs_mat, uint32_t *acc,
                                   int mat_rows, int mat_cols, int bs_mat_cols);

void bitsliced_m_calculate_PS(const uint32_t *P1, const uint32_t *P2,
                              const uint32_t *P3, const unsigned char *S, int m,
                              int v, int o, int k, uint32_t *PS);

void P1_times_O(const mayo_params_t* p, const uint32_t* P1, const unsigned char* O, uint32_t* acc);
void P1P1t_times_O(const mayo_params_t* p, const uint32_t* P1P1t, const unsigned char* O, uint32_t* acc);
void P1_times_Vt(const mayo_params_t* p, const uint32_t* P1, const unsigned char* V, uint32_t* acc);

#endif