
// SPDX-License-Identifier: Apache-2.0

#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include <stdint.h>
#include <mayo.h>
#include <stdint.h>

#if defined(MAYO_AVX) && (M_MAX == 64)
    #include <shuffle_arithmetic_64.h>
#endif
#if defined(MAYO_AVX) && (M_MAX == 96)
    #include <shuffle_arithmetic_96.h>
#endif
#if defined(MAYO_AVX) && (M_MAX == 128)
    #include <shuffle_arithmetic_128.h>
#endif
#if !defined(MAYO_VARIANT) || (defined(MAYO_VARIANT) && (M_MAX == 64))
    #include <arithmetic_64.h>
    #include <arithmetic_96.h>
#endif
#if !defined(MAYO_VARIANT) || (defined(MAYO_VARIANT) && (M_MAX == 96))
    #include <arithmetic_96.h>
    #include <arithmetic_128.h>
#endif
#if !defined(MAYO_VARIANT) || (defined(MAYO_VARIANT) && (M_MAX == 128))
    #include <arithmetic_128.h>
#endif

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

// Calculate P3 = O^T * (P1*O + P2) in KeyGen
void Ot_times_P1O_P2(const mayo_params_t* p, const uint64_t* P1, const unsigned char* O, uint64_t* P1O_P2, uint64_t* P3);

// Calculate Upper in KeyGen
void m_upper(int m_legs, const uint64_t *in, uint64_t *out, int size);

// Calculate acc = (P1+P1^T)*O in expand_sk
void P1P1t_times_O(const mayo_params_t* p, const uint64_t* P1P1t, const unsigned char* O, uint64_t* acc);

// Calculate M=V*L and Y=V*P1*V^T in Sign
void V_times_L__V_times_P1_times_Vt(const mayo_params_t* p, const uint64_t* L, const unsigned char* V, uint64_t* M, const uint64_t* P1, uint64_t* Y);

// Sample solution in Sign
int sample_solution(const mayo_params_t *p, unsigned char *A, const unsigned char *y, const unsigned char *r, unsigned char *x, int k, int o, int m, int A_cols);

// Calculate SPS = S*P*S^T in Verify
void m_calculate_PS_SPS(const uint64_t *P1, const uint64_t *P2, const uint64_t *P3, const unsigned char *S,
                              const int m, const int v, const int o, const int k, uint64_t *SPS);

#endif