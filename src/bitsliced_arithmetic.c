// SPDX-License-Identifier: Apache-2.0

#include "bitsliced_arithmetic.h"
#include <stdalign.h>

// This implements arithmetic for bitsliced vectors of m field elements in Z_2[x]/(x^4+x+1)
// A bitsliced vector is consists of m/32 * 4 consecutive uint32_t's
// the first m/32 uint32_t's are the degree-0 terms, the following m/3 uint32_t's are the degree-1 terms and so on.

void bitsliced_m_vec_mul_add(int m_legs, const uint32_t *in, unsigned char a, uint32_t *acc) {
    const uint32_t *in0 = in;
    const uint32_t *in1 = in + m_legs;
    const uint32_t *in2 = in + 2 * m_legs;
    const uint32_t *in3 = in + 3 * m_legs;

    uint32_t *acc0 = acc;
    uint32_t *acc1 = acc + m_legs;
    uint32_t *acc2 = acc + 2 * m_legs;
    uint32_t *acc3 = acc + 3 * m_legs;

    // terms of a
    uint32_t aa = a;
    uint32_t a0 = -(aa & 1);
    uint32_t a1 = -((aa >> 1) & 1);
    uint32_t a2 = -((aa >> 2) & 1);
    uint32_t a3 = -((aa >> 3) & 1);

    uint32_t x, y, z;
    for (int i = 0; i < m_legs; i++) {
        // deg 0 term of a;
        acc0[i] ^= a0 & in0[i];
        acc1[i] ^= a0 & in1[i];
        acc2[i] ^= a0 & in2[i];
        acc3[i] ^= a0 & in3[i];

        // deg 1 term of a;
        x = in0[i] ^ in3[i];
        acc0[i] ^= a1 & in3[i];
        acc1[i] ^= a1 & x;
        acc2[i] ^= a1 & in1[i];
        acc3[i] ^= a1 & in2[i];

        // deg 2 term of a
        y = in3[i] ^ in2[i];
        acc0[i] ^= a2 & in2[i];
        acc1[i] ^= a2 & y;
        acc2[i] ^= a2 & x;
        acc3[i] ^= a2 & in1[i];

        // deg 3 term of a
        z = in2[i] ^ in1[i];
        acc0[i] ^= a3 & in1[i];
        acc1[i] ^= a3 & z;
        acc2[i] ^= a3 & y;
        acc3[i] ^= a3 & x;
    }
}

inline
static void bitsliced_m_vec_mul_add_x(int m_legs, const uint32_t *in, uint32_t *acc) {
    const uint32_t *in0 = in;
    const uint32_t *in1 = in + m_legs;
    const uint32_t *in2 = in + 2 * m_legs;
    const uint32_t *in3 = in + 3 * m_legs;

    uint32_t *acc0 = acc;
    uint32_t *acc1 = acc + m_legs;
    uint32_t *acc2 = acc + 2 * m_legs;
    uint32_t *acc3 = acc + 3 * m_legs;

    uint32_t x;
    for (int i = 0; i < m_legs; i++) {
        // deg 1 term of a;
        x = in0[i] ^ in3[i];
        acc0[i] ^= in3[i];
        acc1[i] ^= x;
        acc2[i] ^= in1[i];
        acc3[i] ^= in2[i];
    }
}

void bitsliced_m_upper(int m_legs, const uint32_t *in, uint32_t *out, int size) {
    int m_vecs_stored = 0;
    for (int r = 0; r < size; r++) {
        for (int c = r; c < size; c++) {
            bitsliced_m_vec_copy(m_legs, in + m_legs * 4 * (r * size + c), out + m_legs * 4 * m_vecs_stored );
            if (r != c) {
                bitsliced_m_vec_add(m_legs, in + m_legs * 4 * (c * size + r), out + m_legs * 4 * m_vecs_stored );
            }
            m_vecs_stored ++;
        }
    }
}

void bitslice_m_vec(int m_legs, const unsigned char *in, uint32_t *out) {
    uint32_t *out0 = out;
    uint32_t *out1 = out + m_legs;
    uint32_t *out2 = out + 2 * m_legs;
    uint32_t *out3 = out + 3 * m_legs;

    for (int leg = 0; leg < m_legs; leg++) {
        uint32_t d0 = 0, d1 = 0, d2 = 0, d3 = 0;
        for (int i = 31; i >= 0; i--) {
            d0 = (d0 << 1) ^  (in[32 * leg + i] & 1);
            d1 = (d1 << 1) ^ ((in[32 * leg + i] & 2) >> 1);
            d2 = (d2 << 1) ^ ((in[32 * leg + i] & 4) >> 2);
            d3 = (d3 << 1) ^ ((in[32 * leg + i] & 8) >> 3);
        }
        out0[leg] = d0;
        out1[leg] = d1;
        out2[leg] = d2;
        out3[leg] = d3;
    }
}

void unbitslice_m_vec(int m_legs, const uint32_t *in, unsigned char *out) {
    const uint32_t *in0 = in;
    const uint32_t *in1 = in + m_legs;
    const uint32_t *in2 = in + 2 * m_legs;
    const uint32_t *in3 = in + 3 * m_legs;

    for (int leg = 0; leg < m_legs; leg ++) {
        for (int i = 0; i < 32; i++) {
            out[leg * 32 + i] = ((in0[leg] >> i) & 1) ^
                                (((in1[leg] >> i) & 1) << 1) ^
                                (((in2[leg] >> i) & 1) << 2) ^
                                (((in3[leg] >> i) & 1) << 3) ;
        }
    }
}

void P1_times_O(const mayo_params_t* p, const uint32_t* P1, const unsigned char* O, uint32_t* acc){
    #if (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 64 && (O_MAX % 2 == 0))
        mayo_12_P1_times_O(P1,O,acc);
    #elif (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 96 && (O_MAX % 10 == 0))
        mayo_3_P1_times_O(P1,O,acc);
    #elif (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 128 && (O_MAX % 2 == 0))
        mayo_5_P1_times_O(P1,O,acc);
    #else 
        mul_add_bitsliced_m_upper_triangular_mat_x_mat(p->m/32, P1, O, acc, p->n - p->o, p->n - p->o, p->o, 1);
    #endif
}

void P1P1t_times_O(const mayo_params_t* p, const uint32_t* P1P1t, const unsigned char* O, uint32_t* acc){
    #if (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 64 && (O_MAX % 2 == 0))
        mayo_12_P1P1t_times_O(P1P1t,O,acc);
    #elif (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 96 && (O_MAX % 10 == 0))
        mayo_3_P1P1t_times_O(P1P1t,O,acc);
    #elif (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 128 && (O_MAX % 2 == 0))
        mayo_5_P1P1t_times_O(P1P1t,O,acc);
    #else
        mul_add_bitsliced_m_upper_triangular_mat_x_mat(p->m/32, P1P1t, O, acc, p->n - p->o, p->n - p->o, p->o, 0);
    #endif
}

void P1_times_Vt(const mayo_params_t* p, const uint32_t* P1, const unsigned char* V, uint32_t* acc){
    #if (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 64 && (K_MAX == 4 || K_MAX == 9))
        mayo_12_P1_times_Vt(P1,V,acc);
    #elif (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 96 && (K_MAX % 11 == 0))
        mayo_3_P1_times_Vt(P1,V,acc);
    #elif (MAYO_AVX && defined(MAYO_VARIANT) && M_MAX == 128 && (K_MAX % 2 == 0))
        mayo_5_P1_times_Vt(P1,V,acc);
    #else
        mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(p->m/32, P1, V, acc, p->n - p->o, p->n - p->o, p->k, 1);
    #endif
}

// multiplies m bitsliced (possibly upper triangular) matrices with a single matrix and adds result to acc
void mul_add_bitsliced_m_upper_triangular_mat_x_mat(int m_legs, const uint32_t *bs_mat, const unsigned char *mat, uint32_t *acc, int bs_mat_rows, int bs_mat_cols, int mat_cols, int triangular) {

    int bs_mat_entries_used = 0;
    for (int r = 0; r < bs_mat_rows; r++) {
        for (int c = triangular * r; c < bs_mat_cols; c++) {
            for (int k = 0; k < mat_cols; k += 1) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
                bitsliced_64_vec_mul_add((uint64_t *) bs_mat + 4 * bs_mat_entries_used, mat[c * mat_cols + k], (uint64_t *) acc + 4 * (r * mat_cols + k));
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
                bitsliced_96_vec_mul_add(bs_mat + 12 * bs_mat_entries_used, mat[c * mat_cols + k], acc + 12 * (r * mat_cols + k));
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
                bitsliced_128_vec_mul_add((uint64_t *) bs_mat + 8 * bs_mat_entries_used, mat[c * mat_cols + k], (uint64_t *) acc + 8 * (r * mat_cols + k));
#else
                bitsliced_m_vec_mul_add(m_legs, bs_mat + m_legs * 4 * bs_mat_entries_used, mat[c * mat_cols + k], acc + m_legs * 4 * (r * mat_cols + k));
#endif
            }
            bs_mat_entries_used += 1;
        }
    }
}

// multiplies m bitsliced (possibly upper triangular) matrices with the transpose of a single matrix and adds result to acc
void mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(int m_legs, const uint32_t *bs_mat, const unsigned char *mat, uint32_t *acc, int bs_mat_rows, int bs_mat_cols, int mat_rows, int triangular) {

    int bs_mat_entries_used = 0;
    for (int r = 0; r < bs_mat_rows; r++) {
        for (int c = triangular * r; c < bs_mat_cols; c++) {
            for (int k = 0; k < mat_rows; k += 1) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
                bitsliced_64_vec_mul_add((uint64_t *) bs_mat + 4 * bs_mat_entries_used, mat[k * bs_mat_cols + c], (uint64_t *) acc + 4 * (r * mat_rows + k));
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
                bitsliced_96_vec_mul_add(bs_mat + 12 * bs_mat_entries_used, mat[k * bs_mat_cols + c], acc + 12 * (r * mat_rows + k));
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
                bitsliced_128_vec_mul_add((uint64_t *) bs_mat + 8 * bs_mat_entries_used, mat[k * bs_mat_cols + c], (uint64_t *) acc + 8 * (r * mat_rows + k));
#else
                bitsliced_m_vec_mul_add(m_legs, bs_mat + m_legs * 4 * bs_mat_entries_used, mat[k * bs_mat_cols + c], acc + m_legs * 4 * (r * mat_rows + k));
#endif
            }
            bs_mat_entries_used += 1;
        }
    }
}

// multiplies the transpose of a single matrix with m bitsliced matrices and adds result to acc
void mul_add_mat_trans_x_bitsliced_m_mat(int m_legs, const unsigned char *mat, const uint32_t *bs_mat, uint32_t *acc, int mat_rows, int mat_cols, int bs_mat_cols) {

    for (int r = 0; r < mat_cols; r++) {
        for (int c = 0; c < mat_rows; c++) {
            for (int k = 0; k < bs_mat_cols; k += 1) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
                bitsliced_64_vec_mul_add((uint64_t *)bs_mat + 4 * (c * bs_mat_cols + k), mat[c * mat_cols + r], (uint64_t *) acc + 4 * (r * bs_mat_cols + k));
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
                bitsliced_96_vec_mul_add(bs_mat + 12 * (c * bs_mat_cols + k), mat[c * mat_cols + r], acc + 12 * (r * bs_mat_cols + k));
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
                bitsliced_128_vec_mul_add((uint64_t *)bs_mat + 8 * (c * bs_mat_cols + k), mat[c * mat_cols + r], (uint64_t *) acc + 8 * (r * bs_mat_cols + k));
#else
                bitsliced_m_vec_mul_add(m_legs, bs_mat + m_legs * 4 * (c * bs_mat_cols + k), mat[c * mat_cols + r], acc + m_legs * 4 * (r * bs_mat_cols + k));
#endif
            }
        }
    }
}

// multiplies a single matrix with m bitsliced matrices and adds result to acc
void mul_add_mat_x_bitsliced_m_mat(int m_legs, const unsigned char *mat, const uint32_t *bs_mat, uint32_t *acc, int mat_rows, int mat_cols, int bs_mat_cols) {

    for (int r = 0; r < mat_rows; r++) {
        for (int c = 0; c < mat_cols; c++) {
            for (int k = 0; k < bs_mat_cols; k += 1) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
                bitsliced_64_vec_mul_add((uint64_t *)bs_mat + 4 * (c * bs_mat_cols + k), mat[r * mat_cols + c], (uint64_t *) acc + 4 * (r * bs_mat_cols + k));
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
                bitsliced_96_vec_mul_add(bs_mat + 12 * (c * bs_mat_cols + k), mat[r * mat_cols + c], acc + 12 * (r * bs_mat_cols + k));
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
                bitsliced_128_vec_mul_add((uint64_t *)bs_mat + 8 * (c * bs_mat_cols + k), mat[r * mat_cols + c], (uint64_t *) acc + 8 * (r * bs_mat_cols + k));
#else
                bitsliced_m_vec_mul_add(m_legs, bs_mat + m_legs * 4 * (c * bs_mat_cols + k), mat[r * mat_cols + c], acc + m_legs * 4 * (r * bs_mat_cols + k));
#endif
            }
        }
    }
}

void bitsliced_m_multiply_bins(const int m_legs, uint32_t *bins, uint32_t *out) {

    bitsliced_m_vec_add(m_legs, bins + 15 * m_legs * 4, bins + 12 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 15 * m_legs * 4, bins +  3 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 14 * m_legs * 4, bins +  8 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 14 * m_legs * 4, bins +  6 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 13 * m_legs * 4, bins + 10 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 13 * m_legs * 4, bins +  7 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 12 * m_legs * 4, bins +  8 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 12 * m_legs * 4, bins +  4 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 11 * m_legs * 4, bins +  9 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 11 * m_legs * 4, bins +  2 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 10 * m_legs * 4, bins +  8 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 10 * m_legs * 4, bins +  2 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 9 * m_legs * 4, bins +  8 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 9 * m_legs * 4, bins +  1 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 7 * m_legs * 4, bins +  4 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 7 * m_legs * 4, bins +  3 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 6 * m_legs * 4, bins +  4 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 6 * m_legs * 4, bins +  2 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 5 * m_legs * 4, bins +  4 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 5 * m_legs * 4, bins +  1 * m_legs * 4);

    bitsliced_m_vec_add(m_legs, bins + 3 * m_legs * 4, bins +  2 * m_legs * 4);
    bitsliced_m_vec_add(m_legs, bins + 3 * m_legs * 4, bins +  1 * m_legs * 4);

    bitsliced_m_vec_mul_add_x(m_legs, bins + 8 * m_legs * 4, bins + 4 * m_legs * 4);
    bitsliced_m_vec_mul_add_x(m_legs, bins + 4 * m_legs * 4, bins + 2 * m_legs * 4);
    bitsliced_m_vec_mul_add_x(m_legs, bins + 2 * m_legs * 4, bins + 1 * m_legs * 4);

    bitsliced_m_vec_copy(m_legs, bins + 1 * m_legs * 4, out);
}

// compute P * S^t = [ P1  P2 ] * [S1] = [P1*S1 + P2*S2] in bitsliced form
//                   [  0  P3 ]   [S2]   [        P3*S2]
void bitsliced_m_calculate_PS(const uint32_t *P1, const uint32_t *P2, const uint32_t *P3, const unsigned char *S,
                              const int m, const int v, const int o, const int k, uint32_t *PS) {

    const int n = o + v;
#if defined(MAYO_VARIANT) && ((M_MAX == 64) || (M_MAX == 96) || M_MAX == 128)
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

    mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(m_legs, P1, S1, PS, v, v, k, 1); // P1 * S1
    mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(m_legs, P2, S2, PS, v, o, k, 0); // P2 * S2
    mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(m_legs, P3, S2, PS + v*k*m_legs*4, o, o, k, 1); // P3 * S2
    */

    alignas (32) uint32_t accumulator[16 * M_MAX * K_MAX * N_MAX / 8] = {0};
    int P1_used = 0;
    for (int row = 0; row < v; row++) {
        for (int j = row; j < v; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 9)
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 8 ));
#elif defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 4)
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P1 + P1_used * 8), (uint64_t *) (accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 8 ));
#elif defined(MAYO_VARIANT) && (M_MAX == 96) && (K_MAX == 11)
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 9) * 16 + S[9 * n + j] ) * 12 );
            bitsliced_96_vec_add(P1 + P1_used * 12, accumulator + ( (row * k + 10) * 16 + S[10 * n + j] ) * 12 );
#elif defined(MAYO_VARIANT) && (M_MAX == 128) && (K_MAX == 12)
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 9) * 16 + S[9 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 10) * 16 + S[10 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P1 + P1_used * 16), (uint64_t *) (accumulator + ( (row * k + 11) * 16 + S[11 * n + j] ) * 16 ));
#else
            for (int col = 0; col < k; col++) {
                bitsliced_m_vec_add(m_legs, P1 + P1_used * m_legs * 4, accumulator + ( (row * k + col) * 16 + S[col * n + j] )*m_legs * 4 );
            }
#endif
            P1_used ++;
        }

        for (int j = 0; j < o; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 9)
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 0) * 16 + S[(0 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 1) * 16 + S[(1 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 2) * 16 + S[(2 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 3) * 16 + S[(3 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 4) * 16 + S[(4 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 5) * 16 + S[(5 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 6) * 16 + S[(6 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 7) * 16 + S[(7 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 8) * 16 + S[(8 * n) + j + v] ) * 8 ));
#elif defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 4)
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 0) * 16 + S[(0 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 1) * 16 + S[(1 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 2) * 16 + S[(2 * n) + j + v] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P2 + (row * o + j) * 8), (uint64_t *) ( accumulator + ( (row * k + 3) * 16 + S[(3 * n) + j + v] ) * 8 ));
#elif defined(MAYO_VARIANT) && (M_MAX == 96) && (K_MAX == 11)
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 0) * 16 + S[(0 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 1) * 16 + S[(1 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 2) * 16 + S[(2 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 3) * 16 + S[(3 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 4) * 16 + S[(4 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 5) * 16 + S[(5 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 6) * 16 + S[(6 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 7) * 16 + S[(7 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 8) * 16 + S[(8 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 9) * 16 + S[(9 * n) + j + v] ) * 12 );
            bitsliced_96_vec_add(P2 + (row * o + j) * 12, accumulator + ( (row * k + 10) * 16 + S[(10 * n) + j + v] ) * 12 );
#elif defined(MAYO_VARIANT) && (M_MAX == 128) && (K_MAX == 12)
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 0) * 16 + S[(0 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 1) * 16 + S[(1 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 2) * 16 + S[(2 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 3) * 16 + S[(3 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 4) * 16 + S[(4 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 5) * 16 + S[(5 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 6) * 16 + S[(6 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 7) * 16 + S[(7 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 8) * 16 + S[(8 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 9) * 16 + S[(9 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 10) * 16 + S[(10 * n) + j + v] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P2 + (row * o + j) * 16), (uint64_t *) ( accumulator + ( (row * k + 11) * 16 + S[(11 * n) + j + v] ) * 16 ));
#else
            for (int col = 0; col < k; col++) {
                bitsliced_m_vec_add(m_legs, P2 + (row * o + j)*m_legs * 4, accumulator + ( (row * k + col) * 16 + S[(col * n) + j + v] )*m_legs * 4 );
            }
#endif
        }
    }

    int P3_used = 0;
    for (int row = v; row < n; row++) {
        for (int j = row; j < n; j++) {
#if defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 9)
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 8 ));
#elif defined(MAYO_VARIANT) && (M_MAX == 64) && (K_MAX == 4)
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 8 ));
            bitsliced_64_vec_add((uint64_t *) (P3 + P3_used * 8), (uint64_t *) (accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 8 ));
#elif defined(MAYO_VARIANT) && (M_MAX == 96) && (K_MAX == 11)
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 9) * 16 + S[9 * n + j] ) * 12 );
            bitsliced_96_vec_add(P3 + P3_used * 12, accumulator + ( (row * k + 10) * 16 + S[10 * n + j] ) * 12 );
#elif defined(MAYO_VARIANT) && (M_MAX == 128) && (K_MAX == 12)
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 0) * 16 + S[0 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 1) * 16 + S[1 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 2) * 16 + S[2 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 3) * 16 + S[3 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 4) * 16 + S[4 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 5) * 16 + S[5 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 6) * 16 + S[6 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 7) * 16 + S[7 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 8) * 16 + S[8 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 9) * 16 + S[9 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 10) * 16 + S[10 * n + j] ) * 16 ));
            bitsliced_128_vec_add((uint64_t *) (P3 + P3_used * 16), (uint64_t *) (accumulator + ( (row * k + 11) * 16 + S[11 * n + j] ) * 16 ));
#else
            for (int col = 0; col < k; col++) {
                bitsliced_m_vec_add(m_legs, P3 + P3_used * m_legs * 4, accumulator + ( (row * k + col) * 16 + S[col * n + j] )*m_legs * 4 );
            }
#endif
            P3_used ++;
        }
    }

    // multiply stuff according to the bins of the accumulator and add to PS.
    int i = 0;
    while (i < n * k) {
#if defined(MAYO_VARIANT) && (M_MAX == 64)
        bitsliced_64_multiply_bins(accumulator + i * 16 * 8, PS + i * 8);
        i++;
#elif defined(MAYO_VARIANT) && (M_MAX == 96)
        bitsliced_96_multiply_bins(accumulator + i * 16 * 12, PS + i * 12);
        i++;
#elif defined(MAYO_VARIANT) && (M_MAX == 128)
        bitsliced_128_multiply_bins(accumulator + i * 16 * 16, PS + i * 16);
        i++;
#else
        bitsliced_m_multiply_bins(m_legs, accumulator + i * 16 * m_legs * 4, PS + i * m_legs * 4);
        i++;
#endif
    }
}