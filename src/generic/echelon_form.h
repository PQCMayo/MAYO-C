
// SPDX-License-Identifier: Apache-2.0

#ifndef ECHELON_FORM_H
#define ECHELON_FORM_H

#include <stdalign.h>
#include <stdint.h>
#include <mem.h>

#define MAYO_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MAYO_MIN(x, y) (((x) < (y)) ? (x) : (y))

static inline unsigned char
bitsliced_m_extract_element(int m_legs, const uint32_t *in, int index) {
    const uint32_t *in0 = in;
    const uint32_t *in1 = in + m_legs;
    const uint32_t *in2 = in + 2 * m_legs;
    const uint32_t *in3 = in + 3 * m_legs;
    const int leg = index / 32;
    // element order is 0, 8, 16, 24, 1, 9, 17, 25, 2, 10, 18, 26, ...
    const int shift = ((index & 0x7) << 2) + ((index & 0x18) >> 3);

    return ((in0[leg] >> shift) & 1) ^ (((in1[leg] >> shift) & 1) << 1) ^
           (((in2[leg] >> shift) & 1) << 2) ^ (((in3[leg] >> shift) & 1) << 3);
}

static inline void
ef_bitslice_m_vec(int m_legs, const unsigned char *in, uint32_t *out) {
    uint32_t *out0 = out;
    uint32_t *out1 = out + m_legs;
    uint32_t *out2 = out + 2 * m_legs;
    uint32_t *out3 = out + 3 * m_legs;


    uint32_t in32[4];
    for (int leg = 0; leg < m_legs; leg++) {
        for(int i = 0; i < 4; i++) {
            in32[i] = (in[32*leg + i*8 + 0] <<  0) ^ (in[32*leg + i*8 + 1] <<  4) ^ (in[32*leg + i*8 + 2] <<  8) ^ (in[32*leg + i*8 + 3] << 12) ^
                      (in[32*leg + i*8 + 4] << 16) ^ (in[32*leg + i*8 + 5] << 20) ^ (in[32*leg + i*8 + 6] << 24) ^ (in[32*leg + i*8 + 7] << 28);
        }
        out0[leg] = ((in32[0] & 0x11111111) >> 0) ^ ((in32[1] & 0x11111111) << 1) ^ ((in32[2] & 0x11111111) << 2) ^ ((in32[3] & 0x11111111) << 3);
        out1[leg] = ((in32[0] & 0x22222222) >> 1) ^ ((in32[1] & 0x22222222) >> 0) ^ ((in32[2] & 0x22222222) << 1) ^ ((in32[3] & 0x22222222) << 2);
        out2[leg] = ((in32[0] & 0x44444444) >> 2) ^ ((in32[1] & 0x44444444) >> 1) ^ ((in32[2] & 0x44444444) >> 0) ^ ((in32[3] & 0x44444444) << 1);
        out3[leg] = ((in32[0] & 0x88888888) >> 3) ^ ((in32[1] & 0x88888888) >> 2) ^ ((in32[2] & 0x88888888) >> 1) ^ ((in32[3] & 0x88888888) >> 0);
    }
}
static inline void
ef_unbitslice_m_vec(int m_legs, const uint32_t *in, unsigned char *out) {
    const uint32_t *in0 = in;
    const uint32_t *in1 = in + m_legs;
    const uint32_t *in2 = in + 2 * m_legs;
    const uint32_t *in3 = in + 3 * m_legs;

    uint32_t out32[4];
    for (int leg = 0; leg < m_legs; leg ++) {
        out32[0] = ((in0[leg] & 0x11111111) >> 0) ^ ((in1[leg] & 0x11111111) << 1) ^ ((in2[leg] & 0x11111111) << 2) ^ ((in3[leg] & 0x11111111) << 3);
        out32[1] = ((in0[leg] & 0x22222222) >> 1) ^ ((in1[leg] & 0x22222222) >> 0) ^ ((in2[leg] & 0x22222222) << 1) ^ ((in3[leg] & 0x22222222) << 2);
        out32[2] = ((in0[leg] & 0x44444444) >> 2) ^ ((in1[leg] & 0x44444444) >> 1) ^ ((in2[leg] & 0x44444444) >> 0) ^ ((in3[leg] & 0x44444444) << 1);
        out32[3] = ((in0[leg] & 0x88888888) >> 3) ^ ((in1[leg] & 0x88888888) >> 2) ^ ((in2[leg] & 0x88888888) >> 1) ^ ((in3[leg] & 0x88888888) >> 0);


        for(int i = 0; i < 8; i++) {
            out[32*leg + 0*8 + i] = (out32[0] >> (i*4)) & 0xF;
            out[32*leg + 1*8 + i] = (out32[1] >> (i*4)) & 0xF;
            out[32*leg + 2*8 + i] = (out32[2] >> (i*4)) & 0xF;
            out[32*leg + 3*8 + i] = (out32[3] >> (i*4)) & 0xF;
        }
    }
}

// put matrix in row echelon form with ones on first nonzero entries *in
// constant time*
static inline void EF(unsigned char *A, int nrows, int ncols) {

    alignas (32) uint32_t _pivot_row[(K_MAX * O_MAX + 1 + 31) / 32 * 4];
    alignas (32) uint32_t _pivot_row2[(K_MAX * O_MAX + 1 + 31) / 32 * 4];
    alignas (32) uint32_t bitsliced_A[((K_MAX * O_MAX + 1 + 31) / 32) * 4 * M_MAX];

    int legs = (ncols + 31) / 32;

    // bitslice the matrix A
    for (int i = 0; i < nrows; i++) {
        ef_bitslice_m_vec(legs, A + i * ncols, bitsliced_A + i * legs * 4);
    }

    // pivot row is secret, pivot col is not

    unsigned char inverse;
    int pivot_row = 0;
    for (int pivot_col = 0; pivot_col < ncols; pivot_col++) {

        int pivot_row_lower_bound = MAYO_MAX(0, pivot_col + nrows - ncols);
        int pivot_row_upper_bound = MAYO_MIN(nrows - 1, pivot_col);
        // the pivot row is guaranteed to be between these lower and upper bounds if
        // A has full rank

        // zero out pivot row
        for (int i = 0; i < legs * 4; i++) {
            _pivot_row[i] = 0;
            _pivot_row2[i] = 0;
        }

        // try to get a pivot row in constant time
        unsigned char pivot = 0;
        uint32_t pivot_is_zero = -1;
        for (int row = pivot_row_lower_bound;
                row <= MAYO_MIN(nrows - 1, pivot_row_upper_bound + 32); row++) {

            uint32_t is_pivot_row = ~ct_compare_32(row, pivot_row);
            uint32_t below_pivot_row = ct_is_greater_than(row, pivot_row);

            for (int j = 0; j < legs * 4; j++) {
                _pivot_row[j] ^= (is_pivot_row | (below_pivot_row & pivot_is_zero)) &
                                 bitsliced_A[row * legs * 4 + j];
            }
            pivot = bitsliced_m_extract_element(legs, _pivot_row, pivot_col);
            pivot_is_zero = ~ct_compare_32((int) pivot, 0);
        }

        // multiply pivot row by inverse of pivot
        inverse = inverse_f(pivot);
#if defined(MAYO_VARIANT) && (((K_MAX * O_MAX + 1 + 31) / 32) == 3)
        bitsliced_96_vec_mul_add(_pivot_row, inverse, _pivot_row2);
#else
        bitsliced_m_vec_mul_add(legs, _pivot_row, inverse, _pivot_row2);
#endif

        // conditionally write pivot row to the correct row, if there is a nonzero
        // pivot
        for (int row = pivot_row_lower_bound; row <= pivot_row_upper_bound; row++) {
            uint32_t do_copy = ~ct_compare_32(row, pivot_row) & ~pivot_is_zero;
            uint32_t do_not_copy = ~do_copy;
            for (int col = 0; col < legs * 4; col++) {
                bitsliced_A[row * legs * 4 + col] =
                    (do_not_copy & bitsliced_A[row * legs * 4 + col]) +
                    (do_copy & _pivot_row2[col]);
            }
        }

        // eliminate entries below pivot
        for (int row = pivot_row_lower_bound; row < nrows; row++) {
            unsigned char below_pivot = (unsigned char) (ct_is_greater_than(row, pivot_row) & ~pivot_is_zero);
            unsigned char elt_to_elim = bitsliced_m_extract_element(
                                            legs, bitsliced_A + row * legs * 4, pivot_col);

#if defined(MAYO_VARIANT) && (((K_MAX * O_MAX + 1 + 31) / 32) == 3)
            bitsliced_96_vec_mul_add(_pivot_row2, below_pivot & elt_to_elim,
                                    bitsliced_A + row * legs * 4);
#elif defined(MAYO_VARIANT) && (((K_MAX * O_MAX + 1 + 31) / 32) == 4)
            bitsliced_128_vec_mul_add((uint64_t *)_pivot_row2, below_pivot & elt_to_elim,
                                    (uint64_t *)(bitsliced_A + row * legs * 4));                                    
#else
            bitsliced_m_vec_mul_add(legs, _pivot_row2, below_pivot & elt_to_elim,
                                    bitsliced_A + row * legs * 4);

#endif
                            
        }

        pivot_row += (-(int32_t)(~pivot_is_zero));
    }

    unsigned char temp[(O_MAX * K_MAX + 1 + 32) * 100];

    // unbitslice the matrix A
    for (int i = 0; i < nrows; i++) {
        ef_unbitslice_m_vec(legs, bitsliced_A + i * legs * 4, temp);
        for (int j = 0; j < ncols; j++) {
            A[i * ncols + j] = temp[j];
        }
    }

    mayo_secure_clear(temp, K_MAX * O_MAX + 1 + 32);
    mayo_secure_clear(_pivot_row, (K_MAX * O_MAX + 1 + 31) / 32 * 4 * 4);
    mayo_secure_clear(_pivot_row2, (K_MAX * O_MAX + 1 + 31) / 32 * 4 * 4);
}

#endif