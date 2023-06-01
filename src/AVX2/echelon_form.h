// SPDX-License-Identifier: Apache-2.0

#include <immintrin.h>
#include <stdint.h>

#define MAYO_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MAYO_MIN(x, y) (((x) < (y)) ? (x) : (y))


//
// generate multiplication table for '4-bit' variable 'b'. Taken from OV paper!
//
static inline __m256i tbl32_gf16_multab( uint8_t b ) {
    static const unsigned char __gf16_mulbase[128] __attribute__((aligned(32))) = {
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
        0x00, 0x02, 0x04, 0x06, 0x08, 0x0a, 0x0c, 0x0e, 0x03, 0x01, 0x07, 0x05, 0x0b, 0x09, 0x0f, 0x0d, 0x00, 0x02, 0x04, 0x06, 0x08, 0x0a, 0x0c, 0x0e, 0x03, 0x01, 0x07, 0x05, 0x0b, 0x09, 0x0f, 0x0d,
        0x00, 0x04, 0x08, 0x0c, 0x03, 0x07, 0x0b, 0x0f, 0x06, 0x02, 0x0e, 0x0a, 0x05, 0x01, 0x0d, 0x09, 0x00, 0x04, 0x08, 0x0c, 0x03, 0x07, 0x0b, 0x0f, 0x06, 0x02, 0x0e, 0x0a, 0x05, 0x01, 0x0d, 0x09,
        0x00, 0x08, 0x03, 0x0b, 0x06, 0x0e, 0x05, 0x0d, 0x0c, 0x04, 0x0f, 0x07, 0x0a, 0x02, 0x09, 0x01, 0x00, 0x08, 0x03, 0x0b, 0x06, 0x0e, 0x05, 0x0d, 0x0c, 0x04, 0x0f, 0x07, 0x0a, 0x02, 0x09, 0x01
    };


    __m256i bx = _mm256_set1_epi16( b & 0xf );
    __m256i b1 = _mm256_srli_epi16( bx, 1 );

    const __m256i tab0 = _mm256_load_si256((__m256i const *) (__gf16_mulbase + 32 * 0));
    const __m256i tab1 = _mm256_load_si256((__m256i const *) (__gf16_mulbase + 32 * 1));
    const __m256i tab2 = _mm256_load_si256((__m256i const *) (__gf16_mulbase + 32 * 2));
    const __m256i tab3 = _mm256_load_si256((__m256i const *) (__gf16_mulbase + 32 * 3));

    __m256i mask_1  = _mm256_set1_epi16(1);
    __m256i mask_4  = _mm256_set1_epi16(4);
    __m256i mask_0  = _mm256_setzero_si256();

    return ( tab0 & _mm256_cmpgt_epi16( bx & mask_1, mask_0) )
           ^ ( tab1 & _mm256_cmpgt_epi16( b1 & mask_1, mask_0) )
           ^ ( tab2 & _mm256_cmpgt_epi16( bx & mask_4, mask_0) )
           ^ ( tab3 & _mm256_cmpgt_epi16( b1 & mask_4, mask_0) );
}

// put matrix in row echelon form with ones on first nonzero entries *in
// constant time*
static inline void EF(unsigned char *A, int nrows, int ncols) {

    #define AVX_REGS_PER_ROW ((K_MAX * O_MAX + 1 + 31) / 32)
    #define MAX_COLS (AVX_REGS_PER_ROW * 32)

    __m256i _pivot_row[AVX_REGS_PER_ROW];
    __m256i _pivot_row2[AVX_REGS_PER_ROW];
    __m256i A_avx[AVX_REGS_PER_ROW* M_MAX];

    unsigned char* pivot_row_bytes = (unsigned char*) _pivot_row;
    unsigned char* A_bytes = (unsigned char*) A_avx;

    // load A in the tail of AVX2 registers
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++)
        {
            A_bytes[i*MAX_COLS + (MAX_COLS - ncols) + j] = A[ i*ncols + j ];
        }
    }

    // pivot row is secret, pivot col is not

    unsigned char inverse;
    int pivot_row = 0;
    for (int pivot_col = MAX_COLS - ncols; pivot_col < MAX_COLS; pivot_col++) {

        int pivot_col_rounded = pivot_col/32;

        int pivot_row_lower_bound = MAYO_MAX(0, pivot_col + nrows - MAX_COLS);
        int pivot_row_upper_bound = MAYO_MIN(nrows - 1, pivot_col - MAX_COLS + ncols);
        // the pivot row is guaranteed to be between these lower and upper bounds if
        // A has full rank

        // zero out pivot row
        for (int i = 0; i < AVX_REGS_PER_ROW; i++) {
            _pivot_row[i] = _mm256_set1_epi8(0);
            _pivot_row2[i] = _mm256_set1_epi8(0);
        }

        // try to get a pivot row in constant time
        unsigned char pivot = 0;
        uint32_t pivot_is_zero = -1;
        for (int row = pivot_row_lower_bound;
                row <= MAYO_MIN(nrows - 1, pivot_row_upper_bound + 32); row++) {
            uint32_t is_pivot_row = ~ct_compare_32(row, pivot_row);
            uint32_t below_pivot_row = ct_is_greater_than(row, pivot_row);

            __m256i mask = _mm256_set1_epi32( is_pivot_row | (below_pivot_row & pivot_is_zero) );
            for (int j = 0; j < AVX_REGS_PER_ROW; j++) {
                _pivot_row[j] ^= mask & A_avx[row * AVX_REGS_PER_ROW + j];
            }
            pivot = pivot_row_bytes[pivot_col];
            pivot_is_zero = ~ct_compare_32((int) pivot, 0);
        }

        // multiply pivot row by inverse of pivot
        inverse = inverse_f(pivot);
        __m256i inverse_multab = tbl32_gf16_multab(inverse);

        for (int j = pivot_col_rounded; j < AVX_REGS_PER_ROW; j++) {
            _pivot_row2[j] ^= _mm256_shuffle_epi8(inverse_multab, _pivot_row[j]);
        }

        // conditionally write pivot row to the correct row, if there is a nonzero pivot
        for (int row = pivot_row_lower_bound; row <= pivot_row_upper_bound; row++) {
            __m256i mask = _mm256_set1_epi32(~ct_compare_32(row, pivot_row) & ~pivot_is_zero);
            for (int col = pivot_col_rounded; col < AVX_REGS_PER_ROW; col++) { 
                A_avx[row*AVX_REGS_PER_ROW + col] = _mm256_blendv_epi8(A_avx[row*AVX_REGS_PER_ROW + col], _pivot_row2[col], mask);
            }
        }

        // eliminate entries below pivot
        for (int row = pivot_row_lower_bound; row < nrows; row++) {
            unsigned char below_pivot =  (unsigned char) (ct_is_greater_than(row, pivot_row) & ~pivot_is_zero);
            unsigned char elt_to_elim = A_bytes[row*AVX_REGS_PER_ROW*32 + pivot_col];

            __m256i multab = tbl32_gf16_multab(below_pivot & elt_to_elim);    
            for (int j = pivot_col_rounded; j < AVX_REGS_PER_ROW; j++) { 
                A_avx[row*AVX_REGS_PER_ROW + j] ^= _mm256_shuffle_epi8(multab, _pivot_row2[j]);
            }               
        }
        pivot_row += (-(int32_t)(~pivot_is_zero));
    }

    // write the matrix A back
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            A[i * ncols + j] = A_bytes[i*AVX_REGS_PER_ROW*32 + (MAX_COLS - ncols) + j];
        }
    }

    mayo_secure_clear(_pivot_row, AVX_REGS_PER_ROW * 32);
    mayo_secure_clear(_pivot_row2, AVX_REGS_PER_ROW * 32);
    mayo_secure_clear(A_avx, AVX_REGS_PER_ROW * 32 * nrows);
}
