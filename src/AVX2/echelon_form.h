// SPDX-License-Identifier: Apache-2.0

#include <immintrin.h>
#include <stdint.h>
#include <simple_arithmetic.h>
#include <arithmetic_common.h>

// constant-time GF(16) inverse a^14, via a shuffle
static inline unsigned char inverse_f_ct(unsigned char a) {
    const __m128i tab = _mm_setr_epi8(0, 1, 9, 14, 13, 11, 7, 6,
                                      15, 2, 12, 5, 10, 4, 3, 8);
    __m128i idx = _mm_set1_epi8((char)(a & 0xf));
    return (unsigned char)(_mm_cvtsi128_si32(_mm_shuffle_epi8(tab, idx)) & 0xf);
}


// out[i] = 0xff if bit i of bits is set, else 0: byte i takes source byte i / 8, tests bit i % 8
static inline void ef_bits_to_bytes(uint64_t bits, unsigned char out[64]) {
    const __m256i b = _mm256_set1_epi64x((long long) bits);
    const __m256i bit = _mm256_set1_epi64x(0x8040201008040201LL);
    const __m256i lo = _mm256_setr_epi64x(0, 0x0101010101010101LL, 0x0202020202020202LL, 0x0303030303030303LL);
    const __m256i hi = _mm256_add_epi8(lo, _mm256_set1_epi8(4));
    _mm256_storeu_si256((__m256i *) out, _mm256_cmpeq_epi8(_mm256_shuffle_epi8(b, lo) & bit, bit));
    _mm256_storeu_si256((__m256i *) (out + 32), _mm256_cmpeq_epi8(_mm256_shuffle_epi8(b, hi) & bit, bit));
}

#if defined(MAYO_GFNI) && defined(GFNI_ARITHMETIC_H)

static const unsigned char ef_lane_idx[64] __attribute__((aligned(64))) = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63
};

/* GFNI EF, written once over gfv (zmm or ymm). Rows are packed two nibbles per byte and
 * right-aligned in EF_R registers, so the last register is full. Each row update is one
 * affine per register, with the multiplier matrices built in one batch per column. */
#define EF_NROWS M_MAX
#define EF_NCOLS (K_MAX * O_MAX + 1)
#define EF_W ((int) sizeof(gfv))
#define EF_P_BYTES ((EF_NCOLS + 1) / 2)
#define EF_R ((EF_P_BYTES + EF_W - 1) / EF_W)
#define EF_STRIDE (EF_R * EF_W)
#define EF_OFF (EF_STRIDE - EF_P_BYTES)   /* leading zero bytes of each row */
#define EF_TABS_N (((EF_NROWS + 31) / 32) * 32)
#define EF_BUF_N  (((EF_TABS_N + 63) / 64) * 64)

/* bit i set iff lane i of the 64 bytes at col is nonzero / iff base + i >= x (unsigned) */
#if defined(MAYO_GFNI_ZMM)
static inline uint64_t ef_nz_bits(const unsigned char *col) {
    return (uint64_t) _mm512_test_epi8_mask(_mm512_loadu_si512((const void *) col),
                                            _mm512_loadu_si512((const void *) col));
}
static inline uint64_t ef_ge_bits(int base, int x) {
    const __m512i rows = _mm512_add_epi8(_mm512_load_si512((const void *) ef_lane_idx), _mm512_set1_epi8((char) base));
    return (uint64_t) _mm512_cmpge_epu8_mask(rows, _mm512_set1_epi8((char) x));
}
#else
static inline uint64_t ef_nz_bits(const unsigned char *col) {
    const __m256i z = _mm256_setzero_si256();
    return ~((uint64_t) (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) col), z))
           | (uint64_t) (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *) (col + 32)), z)) << 32);
}
static inline uint64_t ef_ge_bits(int base, int x) {
    const __m256i b = _mm256_set1_epi8((char) base), xv = _mm256_set1_epi8((char) x);
    const __m256i r0 = _mm256_add_epi8(_mm256_load_si256((const __m256i *) ef_lane_idx), b);
    const __m256i r1 = _mm256_add_epi8(_mm256_load_si256((const __m256i *) (ef_lane_idx + 32)), b);
    return (uint64_t) (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(r0, xv), r0))
         | (uint64_t) (uint32_t) _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(r1, xv), r1)) << 32;
}
#endif

/* one pivot column. band = the register holding pivot_col, a constant at each call
   site, so the register loops have fixed trip counts and skip the registers before it
   (all zero from this column on). */
static inline __attribute__((always_inline))
void ef_gfni_col(gfv *Apk, gfv *pr, unsigned char *ef_col, unsigned char *ef_sel, uint64_t *ef_tabs,
                 int pivot_col, int *pivot_row_p, const int band) {
    const int nrows = EF_NROWS, ncols = EF_NCOLS;
    const unsigned char *Apk_bytes = (const unsigned char *) Apk;
    const int pivot_row = *pivot_row_p;
    const int pivot_row_lower_bound = MAYO_MAX(0, pivot_col + nrows - ncols);
    const int pivot_row_upper_bound = MAYO_MIN(nrows - 1, pivot_col);
    /* the pivot row is guaranteed to be between these lower and upper bounds if A has full rank */
    const int search_end = MAYO_MIN(nrows - 1, pivot_row_upper_bound + 32);
    const int pcb = EF_OFF + pivot_col / 2;
    const int pcs = (pivot_col & 1) * 4;

    /* the search and the multipliers read the same column; read it once */
    const int elim_rows = nrows - pivot_row_lower_bound;
#if defined(__clang__)
    #pragma clang loop vectorize(disable)
#endif
    for (int i = 0; i < elim_rows; i++) {
        ef_col[i] = (unsigned char) ((Apk_bytes[(pivot_row_lower_bound + i) * EF_STRIDE + pcb] >> pcs) & 0xf);
    }

    /* select every candidate (row >= pivot_row) up to and including the first nonzero
       one, as a lowest set bit isolate; f == 0 selects every candidate. Only masks and
       public shift amounts. */
    const int nsearch = search_end - pivot_row_lower_bound + 1;
    const uint64_t cand = ef_ge_bits(pivot_row_lower_bound, pivot_row) & ((~(uint64_t) 0) >> (64 - nsearch));
    const uint64_t nz = ef_nz_bits(ef_col) & cand;
    const uint64_t f = nz & (~nz + 1);
    ef_bits_to_bytes(cand & (f | (f - 1)), ef_sel);

    for (int n = band; n < EF_R; n++) {
        pr[n] = gfv_zero();
    }
    for (int i = 0; i < nsearch; i++) {
        const gfv mv = gfv_set1_8((char) ef_sel[i]);
        for (int n = band; n < EF_R; n++) {
            pr[n] ^= mv & Apk[(pivot_row_lower_bound + i) * EF_R + n];
        }
    }
    /* the pivot is the pivot_col entry of the xor of the selected rows (public address) */
    const unsigned char pivot = (unsigned char) ((((const unsigned char *) pr)[pcb] >> pcs) & 0xf);
    const uint32_t pivot_is_zero = ~ct_compare_32((int) pivot, 0);

    /* multiply pivot row by inverse of pivot */
    const gfv mat_inv = gfv_set1_64(mayo_gfni_tab(inverse_f_ct(pivot)));
    for (int n = band; n < EF_R; n++) {
        pr[n] = gfv_affine(pr[n], mat_inv);
    }

    /* the multipliers are one column of A: build their matrices in one batch; rows past
       elim_rows build junk that is never read */
    mayo_gfni_tabs(ef_col, ef_tabs, (size_t) ((elim_rows + 31) & ~31));

    /* rows up to the upper bound may hold or lie above the pivot row: conditionally write
       the pivot row there (if the pivot is nonzero), and eliminate only below it, by
       masking the matrix (tab(a) & m == tab(a & m)) */
    int row = pivot_row_lower_bound;
    for (; row <= pivot_row_upper_bound; row++) {
        const gfv wr = gfv_set1_64((uint64_t) (int64_t) (int32_t) (~ct_compare_32(row, pivot_row) & ~pivot_is_zero));
        gfv mat_e = gfv_set1_64(ef_tabs[row - pivot_row_lower_bound] & ct_64_is_greater_than(row, pivot_row));
        __asm__("" : "+v"(mat_e)); // keep the matrix in a register, see GFNI_MACC2
        for (int n = band; n < EF_R; n++) {
            Apk[row * EF_R + n] = gfv_sel(wr, pr[n], Apk[row * EF_R + n]) ^ gfv_affine(pr[n], mat_e);
        }
    }
    /* the rows further down are all below the pivot row */
    for (; row < nrows; row++) {
        gfv mat_e = gfv_set1_64(ef_tabs[row - pivot_row_lower_bound]);
        __asm__("" : "+v"(mat_e)); // keep the matrix in a register, see GFNI_MACC2
        for (int n = band; n < EF_R; n++) {
            Apk[row * EF_R + n] ^= gfv_affine(pr[n], mat_e);
        }
    }

    *pivot_row_p = pivot_row + (-(int32_t)(~pivot_is_zero));
}

/* put matrix in row echelon form with ones on first nonzero entries in constant time;
   output matches the byte-format EF */
static inline void EF(unsigned char *A, int _nrows, int _ncols) {
    (void) _nrows;
    (void) _ncols;
    const int nrows = EF_NROWS, ncols = EF_NCOLS;

    gfv Apk[EF_R * EF_NROWS];
    gfv pr[EF_R];                                   // secret: the pivot row
    unsigned char *Apk_bytes = (unsigned char *) Apk;
    unsigned char ef_col[EF_BUF_N] = {0};           // secret: one column of A
    unsigned char ef_sel[64];                       // secret: which rows are selected
    uint64_t ef_tabs[EF_TABS_N];                    // secret: the column's matrices

    // pack the rows two nibbles per byte after EF_OFF zero bytes; the zero padding
    // stays zero through all row operations
    memset(Apk, 0, sizeof(Apk));
    for (int i = 0; i < nrows; i++) {
        unsigned char *dst = Apk_bytes + i * EF_STRIDE + EF_OFF;
        for (int j = 0; j < ncols / 2; j++) {
            dst[j] = (unsigned char) (A[i * ncols + 2 * j] | (A[i * ncols + 2 * j + 1] << 4));
        }
#if EF_NCOLS % 2 == 1
        dst[ncols / 2] = A[i * ncols + ncols - 1];
#endif
    }

    // pivot row is secret, pivot col is not; one loop per band (register of pivot_col)
    int pivot_row = 0, pivot_col = 0;
#define EF_BAND(b) \
    for (; (b) < EF_R && pivot_col < MAYO_MIN(ncols, 2 * (EF_W * ((b) + 1) - EF_OFF)); pivot_col++) { \
        ef_gfni_col(Apk, pr, ef_col, ef_sel, ef_tabs, pivot_col, &pivot_row, (b)); \
    }
    EF_BAND(0)
    EF_BAND(1)
    EF_BAND(2)
#undef EF_BAND

    // write the matrix A back, one nibble per byte
    for (int i = 0; i < nrows; i++) {
        const unsigned char *src = Apk_bytes + i * EF_STRIDE + EF_OFF;
        for (int j = 0; j < ncols / 2; j++) {
            A[i * ncols + 2 * j] = src[j] & 0xf;
            A[i * ncols + 2 * j + 1] = src[j] >> 4;
        }
#if EF_NCOLS % 2 == 1
        A[i * ncols + ncols - 1] = src[ncols / 2] & 0xf;
#endif
    }

    mayo_secure_clear(Apk, sizeof(Apk));
    mayo_secure_clear(pr, sizeof(pr));
    mayo_secure_clear(ef_col, sizeof(ef_col));
    mayo_secure_clear(ef_sel, sizeof(ef_sel));
    mayo_secure_clear(ef_tabs, sizeof(ef_tabs));
}

#else

/* put matrix in row echelon form with ones on first nonzero entries in constant time*/
static inline void EF(unsigned char *A, int _nrows, int _ncols) {

    (void) _nrows;
    (void) _ncols;

    #define nrows M_MAX
    #define ncols (K_MAX * O_MAX + 1)

    #define AVX_REGS_PER_ROW ((K_MAX * O_MAX + 1 + 31) / 32)
    #define MAX_COLS (AVX_REGS_PER_ROW * 32)

    __m256i _pivot_row[AVX_REGS_PER_ROW];
    __m256i A_avx[AVX_REGS_PER_ROW* M_MAX];

    unsigned char* A_bytes = (unsigned char*) A_avx;

    // load A in the tail of AVX2 registers
    for (int i = 0; i < nrows; i++) {
        memcpy(A_bytes + i*MAX_COLS + (MAX_COLS - ncols), A + i*ncols, ncols);
    }

    // pivot row is secret, pivot col is not
    unsigned char inverse;
    unsigned char ef_val[64];    // secret: the pivot column over the search range
    unsigned char ef_sel[64];    // secret: which of those rows is selected
    int pivot_row = 0;
    int pivot_col = MAYO_MAX(MAX_COLS - ncols,0);
    /* six copies of the loop, one per 32-column band of pivot_col, so that each sees
       pivot_col / 32 as the constant EF_PCR and its register loops get fixed trip counts */
    #define EF_BAND(s) MAYO_MAX(AVX_REGS_PER_ROW - (s) - 1, 0)
    #define EF_PCR EF_BAND(5)
    for (; pivot_col < MAX_COLS-160; pivot_col++) {
        #include "echelon_form_loop.h"
    }
    #undef EF_PCR
    #define EF_PCR EF_BAND(4)
    for (; pivot_col < MAX_COLS-128; pivot_col++) {
        #include "echelon_form_loop.h"
    }
    #undef EF_PCR
    #define EF_PCR EF_BAND(3)
    for (; pivot_col < MAX_COLS-96; pivot_col++) {
        #include "echelon_form_loop.h"
    }
    #undef EF_PCR
    #define EF_PCR EF_BAND(2)
    for (; pivot_col < MAX_COLS-64; pivot_col++) {
        #include "echelon_form_loop.h"
    }
    #undef EF_PCR
    #define EF_PCR EF_BAND(1)
    for (; pivot_col < MAX_COLS-32; pivot_col++) {
        #include "echelon_form_loop.h"
    }
    #undef EF_PCR
    #define EF_PCR EF_BAND(0)
    for (; pivot_col < MAX_COLS; pivot_col++) {
        #include "echelon_form_loop.h"
    }
    #undef EF_PCR
    #undef EF_BAND

    // write the matrix A back
    for (int i = 0; i < nrows; i++) {
        memcpy(A + i*ncols, A_bytes + i*MAX_COLS + (MAX_COLS - ncols), ncols);
    }
    mayo_secure_clear(_pivot_row, AVX_REGS_PER_ROW * 32);
    mayo_secure_clear(ef_val, sizeof(ef_val));
    mayo_secure_clear(ef_sel, sizeof(ef_sel));
    mayo_secure_clear(A_avx, AVX_REGS_PER_ROW * 32 * nrows);
}

#endif // MAYO_GFNI

static const unsigned char ef_lt_mask[64] __attribute__((aligned(32))) = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/* Ar = A * r on the byte-per-nibble matrix (the caller clears the last column
 * of A). Vector GF(16) multiply via mask arithmetic, avoiding the scalar mul_f
 * blocker that is ~10x slower on clang >= 15 / gcc >= 13. */
static inline void EF_mat_vec_mul(const unsigned char *A, const unsigned char *r, unsigned char *Ar) {

    #define MV_COLS (K_MAX * O_MAX + 1)
    #define MV_REGS ((MV_COLS + 31) / 32)

    // r, zero-padded (padding also nulls products with bytes read past the row end),
    // as discrete logs: a * r = alog[(log a + log r) mod 15], zeros via the log 0 sentinel
    const __m256i logt = _mm256_load_si256((__m256i const *) mayo_gf16_log);
    const __m256i alog = _mm256_load_si256((__m256i const *) mayo_gf16_alog);
    unsigned char rpad[MV_REGS * 32];
    for (int i = 0; i < MV_COLS; i++) {
        rpad[i] = r[i];
    }
    for (int i = MV_COLS; i < MV_REGS * 32; i++) {
        rpad[i] = 0;
    }
    __m256i logr[MV_REGS];
    for (int j = 0; j < MV_REGS; j++) {
        logr[j] = _mm256_shuffle_epi8(logt, _mm256_loadu_si256((const __m256i *)(rpad + 32 * j)));
    }

    for (int i = 0; i < M_MAX; i++) {
        __m256i acc = _mm256_setzero_si256();
        for (int j = 0; j < MV_REGS; j++) {
            __m256i av = _mm256_loadu_si256((const __m256i *)(A + i * MV_COLS + 32 * j));
            __m256i idx = _mm256_add_epi8(_mm256_shuffle_epi8(logt, av), logr[j]);
            idx = _mm256_min_epu8(idx, _mm256_sub_epi8(idx, _mm256_set1_epi8(15)));
            acc ^= _mm256_shuffle_epi8(alog, idx);
        }
        // horizontal xor of all bytes
        __m128i h = _mm_xor_si128(_mm256_castsi256_si128(acc), _mm256_extracti128_si256(acc, 1));
        h = _mm_xor_si128(h, _mm_srli_si128(h, 8));
        h = _mm_xor_si128(h, _mm_srli_si128(h, 4));
        h = _mm_xor_si128(h, _mm_srli_si128(h, 2));
        h = _mm_xor_si128(h, _mm_srli_si128(h, 1));
        Ar[i] = (unsigned char)(_mm_cvtsi128_si32(h) & 0xf);
    }

    mayo_secure_clear(rpad, sizeof(rpad));
    mayo_secure_clear(logr, sizeof(logr));

    #undef MV_COLS
    #undef MV_REGS
}

/* constant-time back substitution on [A | y] (one nibble per byte) into x.
 * A's coefficient columns are read-only, so transpose them once and eliminate
 * 32 rows per vpshufb; the candidate column may be y itself (col == k*o), read
 * live. Output matches the scalar version. */
// transpose a 16x16 byte block: dst[j*ds + i] = src[i*ss + j]
static inline void ef_transpose_16x16(const unsigned char *src, int ss,
                                      unsigned char *dst, int ds) {
    __m128i r[16], t[16];
    for (int i = 0; i < 16; i++) r[i] = _mm_loadu_si128((const __m128i *)(src + (size_t) i * ss));
    for (int i = 0; i < 8; i++) {
        t[2*i]     = _mm_unpacklo_epi8(r[2*i], r[2*i+1]);
        t[2*i + 1] = _mm_unpackhi_epi8(r[2*i], r[2*i+1]);
    }
    for (int g = 0; g < 4; g++) {
        r[4*g + 0] = _mm_unpacklo_epi16(t[4*g+0], t[4*g+2]);
        r[4*g + 1] = _mm_unpackhi_epi16(t[4*g+0], t[4*g+2]);
        r[4*g + 2] = _mm_unpacklo_epi16(t[4*g+1], t[4*g+3]);
        r[4*g + 3] = _mm_unpackhi_epi16(t[4*g+1], t[4*g+3]);
    }
    for (int h = 0; h < 2; h++) for (int j = 0; j < 4; j++) {
        t[8*h + 2*j]     = _mm_unpacklo_epi32(r[8*h+j], r[8*h+j+4]);
        t[8*h + 2*j + 1] = _mm_unpackhi_epi32(r[8*h+j], r[8*h+j+4]);
    }
    for (int j = 0; j < 8; j++) {
        r[2*j]     = _mm_unpacklo_epi64(t[j], t[j+8]);
        r[2*j + 1] = _mm_unpackhi_epi64(t[j], t[j+8]);
    }
    for (int i = 0; i < 16; i++) _mm_storeu_si128((__m128i *)(dst + (size_t) i * ds), r[i]);
}

static inline void EF_backsub(const unsigned char *A, unsigned char *x) {

    #define BS_ACOLS (K_MAX * O_MAX + 1)
    #define BS_KO (K_MAX * O_MAX)
    #define BS_CHUNKS ((M_MAX + 31) / 32)
    #define BS_STRIDE (BS_CHUNKS * 32)

    unsigned char y[BS_STRIDE];
    unsigned char T[BS_KO * BS_STRIDE];

#if defined(__clang__)
    #pragma clang loop vectorize(disable)
#endif
    for (int i = 0; i < M_MAX; i++) {
        y[i] = A[i * BS_ACOLS + BS_KO];
    }
    for (int i = M_MAX; i < BS_STRIDE; i++) {
        y[i] = 0;
    }
    /* T = A^t over the first BS_KO columns, in blocks: the strided scalar copy this
       replaces dominated back substitution. Bounds are public. */
    for (int col = 0; col < BS_KO; col++) {
        for (int i = M_MAX; i < BS_STRIDE; i++) {
            T[col * BS_STRIDE + i] = 0;
        }
    }
    {
        /* the last block of each axis is anchored at the end, overlapping its
           predecessor, so no scalar tail is left; the overlap is rewritten with the
           same values. Both axes are >= 16 in every parameter set. */
        for (int bc = 0; bc < (BS_KO + 15) / 16; bc++) {
            const int col = MAYO_MIN(bc * 16, BS_KO - 16);
            for (int br = 0; br < (M_MAX + 15) / 16; br++) {
                const int i = MAYO_MIN(br * 16, M_MAX - 16);
                ef_transpose_16x16(A + (size_t) i * BS_ACOLS + col, BS_ACOLS,
                                   T + (size_t) col * BS_STRIDE + i, BS_STRIDE);
            }
        }
    }

    for (int row = M_MAX - 1; row >= 0; row--) {
        unsigned char finished = 0;
        int col_upper_bound = MAYO_MIN(row + (32 / (M_MAX - row)), BS_KO);
        // the first nonzero entry in row r is between r and col_upper_bound with probability at least ~1-q^{-32}

        for (int col = row; col <= col_upper_bound; col++) {
            // entry of the (echelonized) row at this column; live y for the last column
            unsigned char elt = (col < BS_KO) ? A[row * BS_ACOLS + col] : y[row];
            unsigned char correct_column = ct_compare_8(elt, 0) & ~finished;
            unsigned char u = correct_column & y[row];
            x[col] ^= u;

            __m256i tab = tbl32_gf16_multab2(u);
            const unsigned char *src = (col < BS_KO) ? (T + col * BS_STRIDE) : y;

            // y[i] ^= u * column[i] for i < row (public bound)
            int full = row / 32;
            int rem = row % 32;
            for (int c = 0; c < full; c++) {
                __m256i yv = _mm256_loadu_si256((const __m256i *)(y + 32 * c));
                __m256i cv = _mm256_loadu_si256((const __m256i *)(src + 32 * c));
                _mm256_storeu_si256((__m256i *)(y + 32 * c), yv ^ _mm256_shuffle_epi8(tab, cv));
            }
            if (rem) {
                __m256i mask = _mm256_loadu_si256((const __m256i *)(ef_lt_mask + 32 - rem));
                __m256i yv = _mm256_loadu_si256((const __m256i *)(y + 32 * full));
                __m256i cv = _mm256_loadu_si256((const __m256i *)(src + 32 * full));
                _mm256_storeu_si256((__m256i *)(y + 32 * full), yv ^ (_mm256_shuffle_epi8(tab, cv) & mask));
            }

            finished = finished | correct_column;
        }
    }

    mayo_secure_clear(y, sizeof(y));
    mayo_secure_clear(T, sizeof(T));

    #undef BS_ACOLS
    #undef BS_KO
    #undef BS_CHUNKS
    #undef BS_STRIDE
}

