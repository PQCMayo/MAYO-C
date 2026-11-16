// SPDX-License-Identifier: Apache-2.0

const int pivot_col_rounded = EF_PCR; // == pivot_col / 32 in this band

int pivot_row_lower_bound = MAYO_MAX(0, pivot_col + nrows - MAX_COLS);
int pivot_row_upper_bound = MAYO_MIN(nrows - 1, pivot_col - MAX_COLS + ncols);
/* the pivot row is guaranteed to be between these lower and upper bounds if A has full rank*/

/* zero out pivot row */
for (int i = pivot_col_rounded; i < AVX_REGS_PER_ROW; i++) {
    _pivot_row[i] = _mm256_set1_epi8(0);
}

/* try to get a pivot row in constant time: select every candidate up to and including
   the first nonzero one, as a lowest set bit isolate rather than a prefix AND carried
   one row at a time. Addresses and trip counts come from the public row/pivot_col. */
const int ef_last = MAYO_MIN(nrows - 1, pivot_row_upper_bound + 32);
const int ef_n = ef_last - pivot_row_lower_bound + 1;
for (int i = 0; i < ef_n; i++) {
    ef_val[i] = A_bytes[(pivot_row_lower_bound + i) * AVX_REGS_PER_ROW * 32 + pivot_col];
}
// which candidates are nonzero; ef_val[ef_n ..] is stale but dropped by ef_cand
const __m256i ef_zero = _mm256_setzero_si256();
const uint64_t ef_nz =
      (uint64_t) (uint32_t) ~_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) ef_val), ef_zero))
    | (uint64_t) (uint32_t) ~_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) (ef_val + 32)), ef_zero)) << 32;
/* the candidates, rows >= pivot_row among the ef_n searched, are one run of bits: the
   secret offset pivot_row - pivot_row_lower_bound, clamped at 0 without a branch and
   at most ef_n <= 36, only enters as a shift count */
const int64_t ef_d = (int64_t) pivot_row - pivot_row_lower_bound;
const uint64_t ef_cand = ((~(uint64_t) 0) >> (64 - ef_n)) & ((~(uint64_t) 0) << (ef_d & ~(ef_d >> 63)));
const uint64_t ef_first = (ef_nz & ef_cand) & (~(ef_nz & ef_cand) + 1);
const uint64_t ef_pfx = ef_cand & (ef_first | (ef_first - 1));
ef_bits_to_bytes(ef_pfx, ef_sel);
for (int i = 0; i < ef_n; i++) {
    __m256i mask = _mm256_set1_epi8((char) ef_sel[i]);
    for (int j = pivot_col_rounded; j < AVX_REGS_PER_ROW; j++) {
        _pivot_row[j] ^= mask & A_avx[(pivot_row_lower_bound + i) * AVX_REGS_PER_ROW + j];
    }
}
// the pivot is the pivot_col entry of the xor of the selected rows, read at the public pivot_col
unsigned char pivot = ((unsigned char *) _pivot_row)[pivot_col];
uint32_t pivot_is_zero = ~ct_compare_32((int) pivot, 0);

/* multiply pivot row by inverse of pivot */
inverse = inverse_f_ct(pivot);
__m256i inverse_multab = tbl32_gf16_multab2(inverse);

for (int j = pivot_col_rounded; j < AVX_REGS_PER_ROW; j++) {
    _pivot_row[j] = _mm256_shuffle_epi8(inverse_multab, _pivot_row[j]);
}

/* conditionally write pivot row to the correct row, if there is a nonzero pivot,
   and eliminate entries below pivot. Only rows up to the upper bound can receive the
   pivot row or lie above it: pivot_row <= pivot_col - MAX_COLS + ncols by induction. */
int row = pivot_row_lower_bound;
for (; row <= pivot_row_upper_bound; row++) {
    unsigned char below_pivot =  (unsigned char) ct_64_is_greater_than(row, pivot_row);
    unsigned char elt_to_elim = A_bytes[row*AVX_REGS_PER_ROW*32 + pivot_col];

    __m256i multab = tbl32_gf16_multab2(below_pivot & elt_to_elim);
    __m256i mask = _mm256_set1_epi32( (int)(~ct_compare_32(row, pivot_row) & ~pivot_is_zero) );
    for (int col = pivot_col_rounded; col < AVX_REGS_PER_ROW; col++) {
        A_avx[row*AVX_REGS_PER_ROW + col] = _mm256_blendv_epi8(A_avx[row*AVX_REGS_PER_ROW + col], _pivot_row[col], mask) ^
                                                _mm256_shuffle_epi8(multab, _pivot_row[col]);
    }
}
/* the rows further down are all below the pivot row: two at a time, one table build */
for (; row + 1 < nrows; row += 2) {
    unsigned char e0 = A_bytes[row*AVX_REGS_PER_ROW*32 + pivot_col];
    unsigned char e1 = A_bytes[(row + 1)*AVX_REGS_PER_ROW*32 + pivot_col];
    __m256i mt = tbl32_gf16_multab2_x2(e0, e1);
    __m256i mt0 = _mm256_permute4x64_epi64(mt, 0x44);
    __m256i mt1 = _mm256_permute4x64_epi64(mt, 0xee);
    for (int j = pivot_col_rounded; j < AVX_REGS_PER_ROW; j++) {
        A_avx[row*AVX_REGS_PER_ROW + j] ^= _mm256_shuffle_epi8(mt0, _pivot_row[j]);
        A_avx[(row + 1)*AVX_REGS_PER_ROW + j] ^= _mm256_shuffle_epi8(mt1, _pivot_row[j]);
    }
}
if (row < nrows) {
    __m256i mt = tbl32_gf16_multab2(A_bytes[row*AVX_REGS_PER_ROW*32 + pivot_col]);
    for (int j = pivot_col_rounded; j < AVX_REGS_PER_ROW; j++) {
        A_avx[row*AVX_REGS_PER_ROW + j] ^= _mm256_shuffle_epi8(mt, _pivot_row[j]);
    }
}

pivot_row += (-(int32_t)(~pivot_is_zero));
