
// SPDX-License-Identifier: Apache-2.0

#ifndef ECHELON_FORM_H
#define ECHELON_FORM_H

#include <stdalign.h>
#include <stdint.h>
#include <mem.h>
#include <arithmetic.h>

#define MAYO_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MAYO_MIN(x, y) (((x) < (y)) ? (x) : (y))

static inline unsigned char
m_extract_element(int m_legs, const uint32_t *in, int index) {
    (void) m_legs;
    const int leg = index / 8;
    const int offset = index % 8;

    return (in[leg] >> (offset*4)) & 0xF;
}

static inline void
ef_pack_m_vec(const unsigned char *in, uint32_t *out, int ncols) {
    int i;
    unsigned char *out8 = (unsigned char *)out;
    for(i = 0; i+1 < ncols; i += 2){
        out8[i/2]  = (in[i+0] << 0) | (in[i+1] << 4);
    }
    if (ncols % 2 == 1){
        out8[i/2]  = (in[i+0] << 0);
    }
}

static inline void
ef_unpack_m_vec(int m_legs, const uint32_t *in, unsigned char *out) {
    const unsigned char *in8 = (const unsigned char *)in;
    for(int i = 0; i < m_legs * 32; i += 2){
        out[i]   = (in8[i/2]) & 0xF;
        out[i+1] = (in8[i/2] >> 4);
    }
}


// put matrix in row echelon form with ones on first nonzero entries *in
// constant time*
static inline void EF(unsigned char *A, int nrows, int ncols) {

    alignas (32) uint32_t _pivot_row[(K_MAX * O_MAX + 1 + 31) / 32 * 4];
    alignas (32) uint32_t _pivot_row2[(K_MAX * O_MAX + 1 + 31) / 32 * 4];
    alignas (32) uint32_t packed_A[((K_MAX * O_MAX + 1 + 31) / 32) * 4 * M_MAX];

    int legs = (ncols + 31) / 32;

    // nibbleslice the matrix A
    for (int i = 0; i < nrows; i++) {
        ef_pack_m_vec(A + i * ncols, packed_A + i * legs * 4, ncols);
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
                                 packed_A[row * legs * 4 + j];
            }
            pivot = m_extract_element(legs, _pivot_row, pivot_col);
            pivot_is_zero = ~ct_compare_32((int) pivot, 0);
        }

        // multiply pivot row by inverse of pivot
        inverse = inverse_f(pivot);
#if defined(MAYO_VARIANT) && (((K_MAX * O_MAX + 1 + 31) / 32) == 3)
        vec_mul_add_96((uint64_t *)_pivot_row, inverse, (uint64_t *)_pivot_row2);
#else
        m_vec_mul_add(legs, _pivot_row, inverse, _pivot_row2);
#endif

        // conditionally write pivot row to the correct row, if there is a nonzero
        // pivot
        for (int row = pivot_row_lower_bound; row <= pivot_row_upper_bound; row++) {
            uint32_t do_copy = ~ct_compare_32(row, pivot_row) & ~pivot_is_zero;
            uint32_t do_not_copy = ~do_copy;
            for (int col = 0; col < legs * 4; col++) {
                packed_A[row * legs * 4 + col] =
                    (do_not_copy & packed_A[row * legs * 4 + col]) +
                    (do_copy & _pivot_row2[col]);
            }
        }

        // eliminate entries below pivot
        for (int row = pivot_row_lower_bound; row < nrows; row++) {
            unsigned char below_pivot = (row > pivot_row);
            unsigned char elt_to_elim = m_extract_element(
                                            legs, packed_A + row * legs * 4, pivot_col);
#if defined(MAYO_VARIANT) && (((K_MAX * O_MAX + 1 + 31) / 32) == 3)
            vec_mul_add_96((uint64_t *)_pivot_row2, below_pivot * elt_to_elim,
                                    (uint64_t *)(packed_A + row * legs * 4));
#elif defined(MAYO_VARIANT) && (((K_MAX * O_MAX + 1 + 31) / 32) == 4)
            vec_mul_add_128((uint64_t *)_pivot_row2, below_pivot * elt_to_elim,
                                    (uint64_t *)(packed_A + row * legs * 4));                                    
#else
            m_vec_mul_add(legs, _pivot_row2, below_pivot * elt_to_elim,
                                    packed_A + row * legs * 4);

#endif
                            
        }

        pivot_row += (-(int32_t)(~pivot_is_zero));
    }

    unsigned char temp[(O_MAX * K_MAX + 1 + 32) * 100];

    // unbitslice the matrix A
    for (int i = 0; i < nrows; i++) {
        ef_unpack_m_vec(legs, packed_A + i * legs * 4, temp);
        for (int j = 0; j < ncols; j++) {
            A[i * ncols + j] = temp[j];
        }
    }

    mayo_secure_clear(temp, K_MAX * O_MAX + 1 + 32);
    mayo_secure_clear(_pivot_row, (K_MAX * O_MAX + 1 + 31) / 32 * 4 * 4);
    mayo_secure_clear(_pivot_row2, (K_MAX * O_MAX + 1 + 31) / 32 * 4 * 4);
}

#endif