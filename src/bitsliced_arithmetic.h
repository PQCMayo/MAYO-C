
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

#ifdef ENABLE_PARAMS_DYNAMIC
#define PARAM_m(p) (p->m)
#define PARAM_n(p) (p->n)
#define PARAM_o(p) (p->o)
#define PARAM_v(p) (p->n - p->o)
#define PARAM_A_cols(p) (p->k * p->o + 1)
#define PARAM_k(p) (p->k)
#define PARAM_q(p) (p->q)
#define PARAM_m_bytes(p) (p->m_bytes)
#define PARAM_O_bytes(p) (p->O_bytes)
#define PARAM_v_bytes(p) (p->v_bytes)
#define PARAM_r_bytes(p) (p->r_bytes)
#define PARAM_P1_bytes(p) (p->P1_bytes)
#define PARAM_P2_bytes(p) (p->P2_bytes)
#define PARAM_P3_bytes(p) (p->P3_bytes)
#define PARAM_csk_bytes(p) (p->csk_bytes)
#define PARAM_esk_bytes(p) (p->esk_bytes)
#define PARAM_cpk_bytes(p) (p->cpk_bytes)
#define PARAM_epk_bytes(p) (p->epk_bytes)
#define PARAM_sig_bytes(p) (p->sig_bytes)
#define PARAM_f_tail(p) (p->f_tail)
#define PARAM_salt_bytes(p) (p->salt_bytes)
#define PARAM_sk_seed_bytes(p) (p->sk_seed_bytes)
#define PARAM_digest_bytes(p) (p->digest_bytes)
#define PARAM_pk_seed_bytes(p) (p->pk_seed_bytes)
#elif defined(MAYO_VARIANT)
#define PARAM_m(p) PARAM_NAME(m)
#define PARAM_n(p) PARAM_NAME(n)
#define PARAM_o(p) PARAM_NAME(o)
#define PARAM_v(p) PARAM_NAME(v)
#define PARAM_A_cols(p) PARAM_NAME(A_cols)
#define PARAM_k(p) PARAM_NAME(k)
#define PARAM_q(p) PARAM_NAME(q)
#define PARAM_m_bytes(p) PARAM_NAME(m_bytes)
#define PARAM_O_bytes(p) PARAM_NAME(O_bytes)
#define PARAM_v_bytes(p) PARAM_NAME(v_bytes)
#define PARAM_r_bytes(p) PARAM_NAME(r_bytes)
#define PARAM_P1_bytes(p) PARAM_NAME(P1_bytes)
#define PARAM_P2_bytes(p) PARAM_NAME(P2_bytes)
#define PARAM_P3_bytes(p) PARAM_NAME(P3_bytes)
#define PARAM_csk_bytes(p) PARAM_NAME(csk_bytes)
#define PARAM_esk_bytes(p) PARAM_NAME(esk_bytes)
#define PARAM_cpk_bytes(p) PARAM_NAME(cpk_bytes)
#define PARAM_epk_bytes(p) PARAM_NAME(epk_bytes)
#define PARAM_sig_bytes(p) PARAM_NAME(sig_bytes)
static const unsigned char f_tail[] = PARAM_NAME(f_tail);
#define PARAM_salt_bytes(p) PARAM_NAME(salt_bytes)
#define PARAM_sk_seed_bytes(p) PARAM_NAME(sk_seed_bytes)
#define PARAM_digest_bytes(p) PARAM_NAME(digest_bytes)
#define PARAM_pk_seed_bytes(p) PARAM_NAME(pk_seed_bytes)
#define PARAM_f_tail(p) f_tail
#else
#error "Parameter not specified"
#endif

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

// Calculate SPS = S*P*S^T in Verify
void bitsliced_m_calculate_PS_SPS(const uint32_t *P1, const uint32_t *P2, const uint32_t *P3, const unsigned char *S,
                              const int m, const int v, const int o, const int k, uint32_t *SPS);


void P1_times_O(const mayo_params_t* p, const uint32_t* P1, const unsigned char* O, uint32_t* acc);
void P1P1t_times_O(const mayo_params_t* p, const uint32_t* P1, const unsigned char* O, uint32_t* acc);
void P1_times_Vt(const mayo_params_t* p, const uint32_t* P1, const unsigned char* V, uint32_t* acc);

int sample_solution(const mayo_params_t *p, unsigned char *A,
                           const unsigned char *y, const unsigned char *r,
                           unsigned char *x, int k, int o, int m, int A_cols);

#endif