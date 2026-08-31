// SPDX-License-Identifier: Apache-2.0

#ifndef GFNI_ARITHMETIC_H
#define GFNI_ARITHMETIC_H

#include <stdint.h>
#include <mayo.h>
#include <immintrin.h>
#include <mem.h>
#include <string.h>
#include <stdalign.h>
#include <arithmetic_common.h>
#include <arithmetic_fixed.h>

// GF(16) arithmetic using GFNI + AVX-512. Multiplication by a constant is an
// F2-linear map on the four nibble bits, i.e. an 8x8 bit-matrix that
// vgf2p8affineqb applies to 64 packed bytes at once, so m-vecs stay packed
// (one masked zmm, two for m > 128) and a multiply-accumulate is one affine
// plus one xor.

// P1/P2 are read straight from the packed keystream layout (m/2 bytes per m-vec),
// skipping the unpack spread. PINV turns an m-vec index into the packed byte
// address -- the identity when m % 16 == 0. A load reads M_VEC_LIMBS_MAX limbs,
// so its tail nibbles are next-vec pad (harmless).
// the first producer into VKtmp/SPS assigns, so mayo.c need not pre-zero them
#define MAYO_ACC_STORES_FIRST 1
#define MAYO_PACKED_P1P2 1
#if (M_MAX % 16) == 0
#define PINV(base, vec)      ((const uint64_t *)(base) + (size_t)(vec) * M_VEC_LIMBS_MAX)
#else
#define PINV(base, vec)      ((const uint64_t *)((const uint8_t *)(base) + (size_t)(vec) * (size_t)(M_MAX / 2)))
#endif

typedef uint64_t mayo_multab_t;
#define MAYO_MULTAB_ALIGN _Alignof(uint64_t)
#define MAYO_V_MULTABS_N  (K_MAX * V_MAX)
#define MAYO_S2_MULTABS_N (K_MAX * O_MAX)
#define MAYO_O_MULTABS_N (O_MAX * V_MAX)
// select the multab sign orchestration in arithmetic.h
#define MAYO_HAVE_MULTAB_SIGN 1

// route the shuffle-engine builder names (used by echelon_form.h) to GFNI
#define mayo_O_multabs mayo_gfni_O_multabs
#define mayo_V_multabs mayo_gfni_V_multabs
#define mayo_S1_multabs mayo_gfni_S1_multabs
#define mayo_S2_multabs mayo_gfni_S2_multabs

// The tabs are built as a straight elementwise map so the batched builder applies
// directly. V/S1/S2 are therefore stored in their input (row-major) order rather than
// transposed, and their readers walk them with a stride of V_TAB_STRIDE/S2_TAB_STRIDE.
// The m = 64 lane-fill reads them pre-interleaved instead; see MAYO_ILV_GROUP.
#define V_TAB_STRIDE  V_MAX
#define S2_TAB_STRIDE O_MAX

static
inline void mayo_gfni_O_multabs(const unsigned char *O, uint64_t *O_tabs) {
    mayo_gfni_tabs(O, O_tabs, (size_t) V_MAX * O_MAX);
}


static
inline void mayo_gfni_V_multabs(const unsigned char *V, uint64_t *V_tabs) {
    mayo_gfni_tabs(V, V_tabs, (size_t) K_MAX * V_MAX);
}

static
inline void mayo_gfni_S1_multabs(const unsigned char *S1, uint64_t *S1_tabs) {
    mayo_gfni_tabs(S1, S1_tabs, (size_t) K_MAX * V_MAX);
}

static
inline void mayo_gfni_S2_multabs(const unsigned char *S2, uint64_t *S2_tabs) {
    mayo_gfni_tabs(S2, S2_tabs, (size_t) K_MAX * O_MAX);
}

// register geometry of one m-vec
#if M_VEC_LIMBS_MAX <= 8
#define GFNI_R 1
#define GFNI_TAIL_OFF 0
#define GFNI_TAIL_MASK ((__mmask8)((1u << M_VEC_LIMBS_MAX) - 1))
#else
#define GFNI_R 2
#define GFNI_TAIL_OFF 8
#define GFNI_TAIL_MASK ((__mmask8)((1u << (M_VEC_LIMBS_MAX - 8)) - 1))
#endif

static inline void gfni_load_mvec(const uint64_t *p, __m512i in[GFNI_R]) {
#if GFNI_R == 2
    in[0] = _mm512_loadu_si512((const void *) p);
#endif
    in[GFNI_R - 1] = _mm512_loadu_si512((const void *) (p + GFNI_TAIL_OFF));
}

// Loads are full-width even though an m-vec is only M_VEC_LIMBS_MAX limbs: the masked
// store bounds what is written, so the extra lanes are discarded. gcc schedules the
// masked load badly (2x slower on the products); buffers carry MAYO_MVEC_SLACK limbs
// so the wider read stays in bounds.
// acc ^= t, touching exactly M_VEC_LIMBS_MAX limbs
/* at m = 64 a despread pair fills one register exactly: one store, not two masked halves */
#if M_VEC_LIMBS_MAX == 4 && GFNI_R == 1
#define GFNI_PAIR_STORE_FITS 1
static inline void gfni_xor_store_pair(uint64_t *p, const __m512i x[GFNI_R], const __m512i y[GFNI_R]) {
    const __m512i v = _mm512_shuffle_i32x4(x[0], y[0], 0x44);
    _mm512_storeu_si512((void *) p, _mm512_xor_si512(_mm512_loadu_si512((const void *) p), v));
}
#endif
static inline void gfni_xor_store_mvec(uint64_t *p, const __m512i t[GFNI_R]) {
#if GFNI_R == 2
    _mm512_storeu_si512((void *) p, _mm512_xor_si512(_mm512_loadu_si512((const void *) p), t[0]));
#endif
    _mm512_mask_storeu_epi64((void *) (p + GFNI_TAIL_OFF), GFNI_TAIL_MASK,
                             _mm512_xor_si512(_mm512_loadu_si512((const void *) (p + GFNI_TAIL_OFF)),
                                              t[GFNI_R - 1]));
}

// acc = t (no accumulate). Producers that write every element of their output exactly
// once use this, so the target needs no zeroing beforehand; it also drops a load.
static inline void gfni_store_mvec(uint64_t *p, const __m512i t[GFNI_R]) {
#if GFNI_R == 2
    _mm512_storeu_si512((void *) p, t[0]);
#endif
    _mm512_mask_storeu_epi64((void *) (p + GFNI_TAIL_OFF), GFNI_TAIL_MASK, t[GFNI_R - 1]);
}

// 2x2 GF(16) blocks: one affine can apply any F2-linear map to a byte, so with two
// m-vecs nibble-interleaved (byte j = x's nibble j low, y's high) it computes
// (a*x + b*y, c*x + d*y) -- two inputs and two outputs per affine instead of one.
// Byte j of the matrix qword is output bit 7-j, so out-low sits in bytes 4..7.

// slots, not pairs: an odd k leaves a last slot whose second row is zero, so its
// block computes the leftover output in the low half and nothing in the high
#define KP ((K_MAX + 1) / 2)
// registers an interleaved pair spans: one byte per GF(16) element of the two m-vecs
#define ZP ((M_MAX + 63) / 64)
// the consumers hold one accumulator per output and register, so chunk the outputs to
// keep that inside the register file; ZP == 2 needs no chunking
#define KCH ((24 / ZP) < K_MAX ? (24 / ZP) : K_MAX)
#define OCH ((24 / ZP) < O_MAX ? (24 / ZP) : O_MAX)
#define GFNI_2X2_AM 0x0F0F0F0F00000000ULL
#define GFNI_2X2_BM 0x00000000F0F0F0F0ULL
static inline uint64_t gfni_blk2(uint64_t ta, uint64_t tb, uint64_t tc, uint64_t td) {
    return (ta & GFNI_2X2_AM) | ((tb & GFNI_2X2_AM) << 4)
         | ((tc & GFNI_2X2_BM) >> 4) | (td & GFNI_2X2_BM);
}
static inline void gfni_spread2(const __m512i x[GFNI_R], const __m512i y[GFNI_R], __m512i z[ZP]) {
    const __m512i lo = _mm512_set1_epi8(0x0f);
    for (size_t _n = 0; _n < GFNI_R; _n++) {
        __m512i _x = x[_n], _y = y[_n];
        __m512i _x_shifted = _mm512_srli_epi16(_x, 4), _y_shifted = _mm512_slli_epi16(_y, 4);
        // select low nibble from x and high nibbles from y_shifted
        __m512i z0 = _mm512_ternarylogic_epi64(lo, _x, _y_shifted, 0xCA);
        // select low nibble from x_shifted and high nibbles from y
        __m512i z1 = _mm512_ternarylogic_epi64(lo, _x_shifted, _y, 0xCA);

        if (2 * _n + 1 < ZP){
            z[2 * _n] = z0;
            z[2 * _n + 1] = z1;
        } else {
            z[2 * _n] = _mm512_shuffle_i32x4(z0, z1, 0x44);
        }
    }
}

static inline void gfni_despread2(const __m512i z[ZP], __m512i x[GFNI_R], __m512i y[GFNI_R]) {
    //const __m512i idx = _mm512_setr_epi64(0, 8, 1, 9, 2, 10, 3, 11);
    const __m512i lo = _mm512_set1_epi8(0x0f);
    for (size_t _n = 0; _n < GFNI_R; _n++) {

        __m512i z0 = z[2 * _n], z1;
        if (2 * _n + 1 < ZP){
            z1 = z[2 * _n + 1];
        } else {
            z1 = _mm512_shuffle_i32x4(z0, z0, 0xEE);
        }
        
        __m512i z0_shifted = _mm512_srli_epi16(z0, 4), z1_shifted = _mm512_slli_epi16(z1, 4);
        // select low nibble from x and high nibbles from y_shifted
        x[_n] = _mm512_ternarylogic_epi64(lo, z0, z1_shifted, 0xCA);
        // select low nibble from x_shifted and high nibbles from y
        y[_n] = _mm512_ternarylogic_epi64(lo, z0_shifted, z1, 0xCA);
    }
}
// PS keeps the interleaved pairs, so nothing de-spreads until SPS: one pair is two
// zmm, KP pairs to a row
#define MAYO_VKTMP_N (V_MAX * KP * ZP * 8)
/* as KP, pairing along O_MAX */
#define OP ((O_MAX + 1) / 2)
/* GFNI_BLK_TABLE for reduction-major tables, TABS[col * NOUT + out]. The blocks are
   pre-broadcast, so the affine takes its matrix from memory, one broadcast per ZP */
#define GFNI_BLK_TABLE_RM_BC(B, TABS, NOUT, NCOLS, NP)                                       \
    for (size_t _lp = 0; _lp < (size_t)(NP); _lp++) {                                        \
        const int _l1 = (2 * _lp + 1 < (size_t)(NOUT));                                      \
        for (size_t _cp = 0; _cp < ((size_t)(NCOLS) + 1) / 2; _cp++) {                       \
            const int _c1 = (2 * _cp + 1 < (size_t)(NCOLS));                                 \
            const uint64_t _a = (TABS)[(2 * _cp) * (size_t)(NOUT) + 2 * _lp];                \
            const uint64_t _b = _c1 ? (TABS)[(2 * _cp + 1) * (size_t)(NOUT) + 2 * _lp] : 0;  \
            const uint64_t _c = _l1 ? (TABS)[(2 * _cp) * (size_t)(NOUT) + 2 * _lp + 1] : 0;  \
            const uint64_t _d = (_l1 && _c1)                                                 \
                              ? (TABS)[(2 * _cp + 1) * (size_t)(NOUT) + 2 * _lp + 1] : 0;    \
            (B)[_cp * (size_t)(NP) + _lp] =                                                  \
                _mm512_set1_epi64((long long) gfni_blk2(_a, _b, _c, _d));                    \
        }                                                                                    \
    }
/* GFNI_MACC2 with an explicit slot count, against a pre-broadcast table */
#define GFNI_MACC2_BC(t, z, blk, NP)                                                         \
    for (size_t _p = 0; _p < (size_t)(NP); _p++) {                                           \
        for (size_t _q = 0; _q < ZP; _q++) {                                                 \
            (t)[_p][_q] ^= _mm512_gf2p8affine_epi64_epi8((z)[_q], (blk)[_p], 0);             \
        }                                                                                    \
    }
#define PS2X2_PAIR(base, row, pair) ((base) + ((size_t)(row) * KP + (size_t)(pair)) * (ZP * 8))

/* lambda into column 0 of P1 * V^t. Column 0 is the even half of pair 0 here, so
   spread it against zero and xor the pair slot, the same full-width slot the
   producers store. v m-vector XORs, no multiplications. */
#define MAYO_HAVE_VKT_ADD_LAMBDA 1
static inline void m_vkt_add_lambda(uint64_t *VKtmp, const uint64_t *lambda, int v) {
    __m512i zero[GFNI_R];
    for (size_t q = 0; q < GFNI_R; q++) {
        zero[q] = _mm512_setzero_si512();
    }
    for (int r = 0; r < v; r++) {
        uint64_t *slot = PS2X2_PAIR(VKtmp, r, 0);
        __m512i in[GFNI_R], z[ZP];
        gfni_load_mvec(lambda + (size_t) r * M_VEC_LIMBS_MAX, in);
        gfni_spread2(in, zero, z);
        for (size_t w = 0; w < ZP; w++) {
            _mm512_storeu_si512((void *) (slot + 8 * w),
                                _mm512_xor_si512(_mm512_loadu_si512((const void *) (slot + 8 * w)), z[w]));
        }
    }
}
#define GFNI_BLK_TABLE(B, TABS, STRIDE, NCOLS)                                               \
    for (size_t _lp = 0; _lp < KP; _lp++) {                                                  \
        const uint64_t *_t0 = (TABS) + (2 * _lp) * (size_t)(STRIDE);                         \
        const uint64_t *_t1 = (2 * _lp + 1 < K_MAX)                                          \
                            ? (TABS) + (2 * _lp + 1) * (size_t)(STRIDE) : (const uint64_t *) 0; \
        for (size_t _cp = 0; _cp < ((size_t)(NCOLS) + 1) / 2; _cp++) {                       \
            const int _c1 = (2 * _cp + 1 < (size_t)(NCOLS));                                 \
            (B)[_cp * KP + _lp] = gfni_blk2(_t0[2 * _cp], _c1 ? _t0[2 * _cp + 1] : 0,         \
                                            _t1 ? _t1[2 * _cp] : 0,                          \
                                            (_t1 && _c1) ? _t1[2 * _cp + 1] : 0);             \
        }                                                                                    \
    }
// one 2x2 accumulate over both halves of an interleaved pair. The empty asm keeps the
// matrix in a register: clang 14..19 miscompiles the folded m64bcst operand in unrolled
// loops, stepping its displacement by 64 bytes instead of 8, and the broadcast is free.
#define GFNI_MACC2(t, z, blk)                                                                \
    for (size_t _p = 0; _p < KP; _p++) {                                                     \
        __m512i _m = _mm512_set1_epi64((long long) (blk)[_p]);                               \
        __asm__("" : "+v"(_m));                                                              \
        for (size_t _q = 0; _q < ZP; _q++) {                                                 \
            (t)[_p][_q] ^= _mm512_gf2p8affine_epi64_epi8((z)[_q], _m, 0);                    \
        }                                                                                    \
    }
// The pair is interleaved on the free index, and the scalar here depends on the
// contraction row and the output, not on which half -- so both halves take the same
// matrix and this side needs no 2x2. SPS is small, so it is where the de-spread goes.
#define GFNI_SPS_ROW(PSB, ROWS, TABS, STRIDE, CP, L0, LN, T)                                         \
    for (size_t _r = 0; _r < (size_t)(ROWS); _r++) {                                         \
        const uint64_t *_p = PS2X2_PAIR(PSB, _r, CP);                                        \
        __m512i _z[ZP];                                                                      \
        for (size_t _q = 0; _q < ZP; _q++) {                                                 \
            _z[_q] = _mm512_loadu_si512((const void *) (_p + 8 * _q));                       \
        }                                                                                    \
        for (size_t _l = 0; _l < (size_t)(LN); _l++) {                                       \
            __m512i _m = _mm512_set1_epi64((long long)                                       \
                             (TABS)[((size_t)(L0) + _l) * (size_t)(STRIDE) + _r]);           \
            __asm__("" : "+v"(_m));                                                          \
            for (size_t _q = 0; _q < ZP; _q++) {                                             \
                (T)[_l][_q] ^= _mm512_gf2p8affine_epi64_epi8(_z[_q], _m, 0);                 \
            }                                                                                \
        }                                                                                    \
    }




// P1 * O -> P1: v x v (upper triangular), O: v x o
static
inline void P1_times_O(const uint64_t *P1, const uint64_t *O_tabs, uint64_t *acc) {
    __m512i BO[((V_MAX + 1) / 2) * OP];
    GFNI_BLK_TABLE_RM_BC(BO, O_tabs, O_MAX, V_MAX, OP);
    size_t cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++) {
        __m512i t[OP][ZP];
        for (size_t lp = 0; lp < OP; lp++) {
            for (size_t q = 0; q < ZP; q++) t[lp][q] = _mm512_setzero_si512();
        }
        size_t c = r;
        if (c & 1) {            /* an odd start puts the first column in the pair's second half */
            __m512i in[GFNI_R], z[ZP];
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_load_mvec(PINV(P1, cols_used), in);
            cols_used++;
            gfni_spread2(zero, in, z);
            GFNI_MACC2_BC(t, z, BO + (c / 2) * OP, OP);
            c++;
        }
        for (; c + 1 < V_MAX; c += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, cols_used), a);
            gfni_load_mvec(PINV(P1, cols_used + 1), b);
            cols_used += 2;
            gfni_spread2(a, b, z);
            GFNI_MACC2_BC(t, z, BO + (c / 2) * OP, OP);
        }
        if (c < V_MAX) {        /* odd V_MAX */
            __m512i in[GFNI_R], z[ZP];
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_load_mvec(PINV(P1, cols_used), in);
            cols_used++;
            gfni_spread2(in, zero, z);
            GFNI_MACC2_BC(t, z, BO + (c / 2) * OP, OP);
        }
        for (size_t lp = 0; lp < OP; lp++) {
            __m512i x[GFNI_R], y[GFNI_R];
            gfni_despread2(t[lp], x, y);
#ifdef GFNI_PAIR_STORE_FITS
            if (2 * lp + 1 < O_MAX) {
                gfni_xor_store_pair(acc + (r * O_MAX + 2 * lp) * M_VEC_LIMBS_MAX, x, y);
                continue;
            }
#endif
            gfni_xor_store_mvec(acc + (r * O_MAX + 2 * lp) * M_VEC_LIMBS_MAX, x);
            if (2 * lp + 1 < O_MAX) {
                gfni_xor_store_mvec(acc + (r * O_MAX + 2 * lp + 1) * M_VEC_LIMBS_MAX, y);
            }
        }
    }
}

static
inline void Ot_times_P1O_P2(const uint64_t *P1O_P2, const uint64_t *O_tabs, uint64_t *acc) {
    __m512i BO[((V_MAX + 1) / 2) * OP];
    GFNI_BLK_TABLE_RM_BC(BO, O_tabs, O_MAX, V_MAX, OP);
    for (size_t c = 0; c < O_MAX; c++) {
        __m512i t[OP][ZP];
        for (size_t lp = 0; lp < OP; lp++) {
            for (size_t q = 0; q < ZP; q++) t[lp][q] = _mm512_setzero_si512();
        }
        size_t r = 0;
        for (; r + 1 < V_MAX; r += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(P1O_P2 + (r * O_MAX + c) * M_VEC_LIMBS_MAX, a);
            gfni_load_mvec(P1O_P2 + ((r + 1) * O_MAX + c) * M_VEC_LIMBS_MAX, b);
            gfni_spread2(a, b, z);
            GFNI_MACC2_BC(t, z, BO + (r / 2) * OP, OP);
        }
        if (r < V_MAX) {        /* odd V_MAX */
            __m512i in[GFNI_R], z[ZP];
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_load_mvec(P1O_P2 + (r * O_MAX + c) * M_VEC_LIMBS_MAX, in);
            gfni_spread2(in, zero, z);
            GFNI_MACC2_BC(t, z, BO + (r / 2) * OP, OP);
        }
        for (size_t lp = 0; lp < OP; lp++) {
            __m512i x[GFNI_R], y[GFNI_R];
            gfni_despread2(t[lp], x, y);
            gfni_xor_store_mvec(acc + ((2 * lp) * O_MAX + c) * M_VEC_LIMBS_MAX, x);
            if (2 * lp + 1 < O_MAX) {
                gfni_xor_store_mvec(acc + ((2 * lp + 1) * O_MAX + c) * M_VEC_LIMBS_MAX, y);
            }
        }
    }
}

static
inline void P1P1t_times_O(const mayo_params_t *p, const uint64_t *P1, const unsigned char *O, uint64_t *acc) {
    (void) p;
    uint64_t O_tabs[MAYO_O_MULTABS_N];
    mayo_O_multabs(O, O_tabs);
    __m512i BO[((V_MAX + 1) / 2) * OP];
    GFNI_BLK_TABLE_RM_BC(BO, O_tabs, O_MAX, V_MAX, OP);

    __m512i zero[GFNI_R];
    for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();

    size_t diag = 0;            /* packed index of (r, r) */
    for (size_t r = 0; r < V_MAX; r++) {
        __m512i t[OP][ZP];
        for (size_t lp = 0; lp < OP; lp++) {
            for (size_t q = 0; q < ZP; q++) t[lp][q] = _mm512_setzero_si512();
        }
        /* columns left of r come down column r, those right of r along row r; column r
           drops out, P1 + P1t having a zero diagonal. The split keeps its pair out of
           both runs. */
        const size_t rp = r & ~(size_t)1;
        size_t pos = r;         /* packed index of (c, r) */
        size_t c = 0;
        for (; c < rp; c += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, pos), a);
            pos += (V_MAX - c - 1);
            gfni_load_mvec(PINV(P1, pos), b);
            pos += (V_MAX - c - 2);
            gfni_spread2(a, b, z);
            GFNI_MACC2_BC(t, z, BO + (c / 2) * OP, OP);
        }
        {                       /* the pair holding the diagonal: one column, or none */
            __m512i in[GFNI_R], z[ZP];
            if (r & 1) {
                gfni_load_mvec(PINV(P1, pos), in);   /* (r - 1, r) */
                gfni_spread2(in, zero, z);
                GFNI_MACC2_BC(t, z, BO + (rp / 2) * OP, OP);
            } else if (r + 1 < V_MAX) {
                gfni_load_mvec(PINV(P1, diag + 1), in);
                gfni_spread2(zero, in, z);
                GFNI_MACC2_BC(t, z, BO + (rp / 2) * OP, OP);
            }
            c = rp + 2;
        }
        for (; c + 1 < V_MAX; c += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, diag + (c - r)), a);
            gfni_load_mvec(PINV(P1, diag + (c + 1 - r)), b);
            gfni_spread2(a, b, z);
            GFNI_MACC2_BC(t, z, BO + (c / 2) * OP, OP);
        }
        if (c < V_MAX) {        /* odd V_MAX */
            __m512i in[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, diag + (c - r)), in);
            gfni_spread2(in, zero, z);
            GFNI_MACC2_BC(t, z, BO + (c / 2) * OP, OP);
        }
        for (size_t lp = 0; lp < OP; lp++) {
            __m512i x[GFNI_R], y[GFNI_R];
            gfni_despread2(t[lp], x, y);
#ifdef GFNI_PAIR_STORE_FITS
            if (2 * lp + 1 < O_MAX) {
                gfni_xor_store_pair(acc + (r * O_MAX + 2 * lp) * M_VEC_LIMBS_MAX, x, y);
                continue;
            }
#endif
            gfni_xor_store_mvec(acc + (r * O_MAX + 2 * lp) * M_VEC_LIMBS_MAX, x);
            if (2 * lp + 1 < O_MAX) {
                gfni_xor_store_mvec(acc + (r * O_MAX + 2 * lp + 1) * M_VEC_LIMBS_MAX, y);
            }
        }
        diag += (V_MAX - r);
    }
    mayo_secure_clear(O_tabs, sizeof(O_tabs));
}

static
inline void Vt_times_L(const uint64_t *L, const uint64_t *V_tabs, uint64_t *acc) {
    /* L is in the plain layout, so the pair is spread here, not loaded interleaved */
    uint64_t BV[((V_MAX + 1) / 2) * KP];
    GFNI_BLK_TABLE(BV, V_tabs, V_TAB_STRIDE, V_MAX);
    for (size_t c = 0; c < O_MAX; c++) {
        __m512i t[KP][ZP];
        for (size_t lp = 0; lp < KP; lp++) {
            for (size_t q = 0; q < ZP; q++) t[lp][q] = _mm512_setzero_si512();
        }
        size_t r = 0;
        for (; r + 1 < V_MAX; r += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(L + (r * O_MAX + c) * M_VEC_LIMBS_MAX, a);
            gfni_load_mvec(L + ((r + 1) * O_MAX + c) * M_VEC_LIMBS_MAX, b);
            gfni_spread2(a, b, z);
            GFNI_MACC2(t, z, BV + (r / 2) * KP);
        }
        if (r < V_MAX) {            /* odd V_MAX: the last row pairs with zero */
            __m512i in[GFNI_R], z[ZP], zero[GFNI_R];
            gfni_load_mvec(L + (r * O_MAX + c) * M_VEC_LIMBS_MAX, in);
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_spread2(in, zero, z);
            GFNI_MACC2(t, z, BV + (r / 2) * KP);
        }
        for (size_t lp = 0; lp < KP; lp++) {
            __m512i x[GFNI_R], y[GFNI_R];
            gfni_despread2(t[lp], x, y);
            gfni_xor_store_mvec(acc + ((2 * lp) * O_MAX + c) * M_VEC_LIMBS_MAX, x);
            if (2 * lp + 1 < K_MAX) {
                gfni_xor_store_mvec(acc + ((2 * lp + 1) * O_MAX + c) * M_VEC_LIMBS_MAX, y);
            }
        }
    }
}

static
inline void Vt_times_Pv(const uint64_t *Pv, const uint64_t *V_tabs, uint64_t *acc) {
    for (size_t cp = 0; cp < KP; cp++) {
        for (size_t l0 = 0; l0 < K_MAX; l0 += KCH) {
            const size_t ln = (K_MAX - l0 < (size_t) KCH) ? K_MAX - l0 : (size_t) KCH;
            __m512i t[KCH][ZP];
            for (size_t li = 0; li < ln; li++) {
                for (size_t q = 0; q < ZP; q++) t[li][q] = _mm512_setzero_si512();
            }
            GFNI_SPS_ROW(Pv, V_MAX, V_tabs, V_TAB_STRIDE, cp, l0, ln, t);
            for (size_t li = 0; li < ln; li++) {
                __m512i x[GFNI_R], y[GFNI_R];
                gfni_despread2(t[li], x, y);
                gfni_store_mvec(acc + ((l0 + li) * K_MAX + 2 * cp) * M_VEC_LIMBS_MAX, x);
                if (2 * cp + 1 < K_MAX) gfni_store_mvec(acc + ((l0 + li) * K_MAX + 2 * cp + 1) * M_VEC_LIMBS_MAX, y);
            }
        }
    }
}

static
inline void P1_times_Vt(const uint64_t *P1, const uint64_t *V_tabs, uint64_t *acc) {
    uint64_t BV[((V_MAX + 1) / 2) * KP];
    GFNI_BLK_TABLE(BV, V_tabs, V_TAB_STRIDE, V_MAX);
    size_t cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++) {
        __m512i t[KP][ZP];
        for (size_t lp = 0; lp < KP; lp++) {
            for (size_t q = 0; q < ZP; q++) t[lp][q] = _mm512_setzero_si512();
        }
        size_t c = r;
        if (c & 1) {
            __m512i in[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, cols_used), in);
            cols_used++;
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_spread2(zero, in, z);
            GFNI_MACC2(t, z, BV + (c / 2) * KP);
            c++;
        }
        for (; c + 1 < V_MAX; c += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, cols_used), a);
            gfni_load_mvec(PINV(P1, cols_used + 1), b);
            cols_used += 2;
            gfni_spread2(a, b, z);
            GFNI_MACC2(t, z, BV + (c / 2) * KP);
        }
        if (c < V_MAX) {        /* odd V_MAX: a last column with no partner */
            __m512i in[GFNI_R], z[ZP], zero[GFNI_R];
            gfni_load_mvec(PINV(P1, cols_used), in);
            cols_used++;
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_spread2(in, zero, z);
            GFNI_MACC2(t, z, BV + (c / 2) * KP);
        }
        for (size_t lp = 0; lp < KP; lp++) {
            uint64_t *q = PS2X2_PAIR(acc, r, lp);
            for (size_t w = 0; w < ZP; w++) _mm512_storeu_si512((void *) (q + 8 * w), t[lp][w]);
        }
    }
}

// acc += P1^t * V^t, where P1 is stored as an upper-triangular matrix.
static
inline void P1t_times_Vt(const uint64_t *P1, const uint64_t *V_tabs, uint64_t *acc) {
    uint64_t BV[((V_MAX + 1) / 2) * KP];
    GFNI_BLK_TABLE(BV, V_tabs, V_TAB_STRIDE, V_MAX);
    for (size_t r = 0; r < V_MAX; r++) {
        __m512i t[KP][ZP];
        for (size_t lp = 0; lp < KP; lp++) {
            for (size_t q = 0; q < ZP; q++) t[lp][q] = _mm512_setzero_si512();
        }
        /* the walk lands exactly on the diagonal element, so one loop covers c = 0..r */
        size_t pos = r, c = 0;
        for (; c + 1 <= r; c += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, pos), a);
            pos += (V_MAX - c - 1);
            gfni_load_mvec(PINV(P1, pos), b);
            pos += (V_MAX - c - 2);
            gfni_spread2(a, b, z);
            GFNI_MACC2(t, z, BV + (c / 2) * KP);
        }
        if (c <= r) {
            __m512i in[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, pos), in);
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_spread2(in, zero, z);
            GFNI_MACC2(t, z, BV + (c / 2) * KP);
        }
        for (size_t lp = 0; lp < KP; lp++) {
            uint64_t *q = PS2X2_PAIR(acc, r, lp);
            for (size_t w = 0; w < ZP; w++) {
                _mm512_storeu_si512((void *) (q + 8 * w),
                    _mm512_loadu_si512((const void *) (q + 8 * w)) ^ t[lp][w]);
            }
        }
    }
}

/* acc[c] += sum_r O[r][c] * lambda[r]: the one column of U^t * O that the
   materialized-L signing path has to supply itself. Same shape as Ut_times_O
   below with a single column pair and only its even half live, so lambda is
   read unspread and nothing has to be de-spread on the way out. */
#define MAYO_HAVE_LAMBDA_TIMES_O 1
static inline void m_lambda_times_O(const uint64_t *lambda, const mayo_multab_t *O_tabs,
                                    uint64_t *acc) {
    for (size_t l0 = 0; l0 < O_MAX; l0 += OCH) {
        const size_t ln = (O_MAX - l0 < (size_t) OCH) ? O_MAX - l0 : (size_t) OCH;
        __m512i t[OCH][GFNI_R];
        for (size_t li = 0; li < ln; li++) {
            for (size_t q = 0; q < GFNI_R; q++) {
                t[li][q] = _mm512_setzero_si512();
            }
        }
        for (size_t r = 0; r < V_MAX; r++) {
            __m512i z[GFNI_R];
            gfni_load_mvec(lambda + r * M_VEC_LIMBS_MAX, z);
            for (size_t li = 0; li < ln; li++) {
                __m512i m = _mm512_set1_epi64((long long) O_tabs[O_MAX * r + l0 + li]);
                __asm__("" : "+v"(m));
                for (size_t q = 0; q < GFNI_R; q++) {
                    t[li][q] ^= _mm512_gf2p8affine_epi64_epi8(z[q], m, 0);
                }
            }
        }
        for (size_t li = 0; li < ln; li++) {
            gfni_xor_store_mvec(acc + (l0 + li) * M_VEC_LIMBS_MAX, t[li]);
        }
    }
}

// acc += U^t * O, where U is a v-by-k matrix of m-vectors.
static
inline void Ut_times_O(const uint64_t *U, const uint64_t *O_tabs, uint64_t *acc) {
    for (size_t cp = 0; cp < KP; cp++) {
        for (size_t l0 = 0; l0 < O_MAX; l0 += OCH) {
            const size_t ln = (O_MAX - l0 < (size_t) OCH) ? O_MAX - l0 : (size_t) OCH;
            __m512i t[OCH][ZP];
            for (size_t li = 0; li < ln; li++) {
                for (size_t q = 0; q < ZP; q++) t[li][q] = _mm512_setzero_si512();
            }
            for (size_t r = 0; r < V_MAX; r++) {
                const uint64_t *p = PS2X2_PAIR(U, r, cp);
                __m512i z[ZP];
                for (size_t q = 0; q < ZP; q++) z[q] = _mm512_loadu_si512((const void *) (p + 8 * q));
                for (size_t li = 0; li < ln; li++) {
                    __m512i m = _mm512_set1_epi64((long long) O_tabs[O_MAX * r + l0 + li]);
                    __asm__("" : "+v"(m));
                    for (size_t q = 0; q < ZP; q++) {
                        t[li][q] ^= _mm512_gf2p8affine_epi64_epi8(z[q], m, 0);
                    }
                }
            }
            for (size_t li = 0; li < ln; li++) {
                __m512i x[GFNI_R], y[GFNI_R];
                gfni_despread2(t[li], x, y);
                gfni_xor_store_mvec(acc + ((2 * cp) * O_MAX + l0 + li) * M_VEC_LIMBS_MAX, x);
                if (2 * cp + 1 < K_MAX) {
                    gfni_xor_store_mvec(acc + ((2 * cp + 1) * O_MAX + l0 + li) * M_VEC_LIMBS_MAX, y);
                }
            }
        }
    }
}


static
inline void S1t_times_PS1(const uint64_t *_PS1, const uint64_t *S1_tabs, uint64_t *_acc) {
    for (size_t cp = 0; cp < KP; cp++) {
        for (size_t l0 = 0; l0 < K_MAX; l0 += KCH) {
            const size_t ln = (K_MAX - l0 < (size_t) KCH) ? K_MAX - l0 : (size_t) KCH;
            __m512i t[KCH][ZP];
            for (size_t li = 0; li < ln; li++) {
                for (size_t q = 0; q < ZP; q++) t[li][q] = _mm512_setzero_si512();
            }
            GFNI_SPS_ROW(_PS1, V_MAX, S1_tabs, V_TAB_STRIDE, cp, l0, ln, t);
            for (size_t li = 0; li < ln; li++) {
                __m512i x[GFNI_R], y[GFNI_R];
                gfni_despread2(t[li], x, y);
                gfni_store_mvec(_acc + ((l0 + li) * K_MAX + 2 * cp) * M_VEC_LIMBS_MAX, x);
                if (2 * cp + 1 < K_MAX) {
                    gfni_store_mvec(_acc + ((l0 + li) * K_MAX + 2 * cp + 1) * M_VEC_LIMBS_MAX, y);
                }
            }
        }
    }
}

static
inline void S2t_times_PS2(const uint64_t *PS2, const uint64_t *S2_tabs, uint64_t *acc) {
    for (size_t cp = 0; cp < KP; cp++) {
        for (size_t l0 = 0; l0 < K_MAX; l0 += KCH) {
            const size_t ln = (K_MAX - l0 < (size_t) KCH) ? K_MAX - l0 : (size_t) KCH;
            __m512i t[KCH][ZP];
            for (size_t li = 0; li < ln; li++) {
                for (size_t q = 0; q < ZP; q++) t[li][q] = _mm512_setzero_si512();
            }
            GFNI_SPS_ROW(PS2, O_MAX, S2_tabs, S2_TAB_STRIDE, cp, l0, ln, t);
            for (size_t li = 0; li < ln; li++) {
                __m512i x[GFNI_R], y[GFNI_R];
                gfni_despread2(t[li], x, y);
                gfni_xor_store_mvec(acc + ((l0 + li) * K_MAX + 2 * cp) * M_VEC_LIMBS_MAX, x);
                if (2 * cp + 1 < K_MAX) gfni_xor_store_mvec(acc + ((l0 + li) * K_MAX + 2 * cp + 1) * M_VEC_LIMBS_MAX, y);
            }
        }
    }
}

// P2*S2 -> P2: v x o, S2: o x k
static
inline void P1_times_S1_plus_P2_times_S2(const uint64_t *P1, const uint64_t *P2, const uint64_t *S1_tabs,
                                         const uint64_t *S2_tabs, uint64_t *acc) {
    /* One affine per input pair and output pair. The blocks come straight from the
       scalar matrices; an odd start spreads (0, x) so the pair entry still applies. */
    uint64_t B1[((V_MAX + 1) / 2) * KP], B2[((O_MAX + 1) / 2) * KP];
    GFNI_BLK_TABLE(B1, S1_tabs, V_TAB_STRIDE, V_MAX);
    GFNI_BLK_TABLE(B2, S2_tabs, S2_TAB_STRIDE, O_MAX);

    size_t P1_vec = 0;
    for (size_t r = 0; r < V_MAX; r++) {
        __m512i t[KP][ZP];
        for (size_t lp = 0; lp < KP; lp++) {
            for (size_t q = 0; q < ZP; q++) t[lp][q] = _mm512_setzero_si512();
        }

        // P1 * S1
        size_t c = r;
        if (c & 1) {
            __m512i in[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, P1_vec), in);
            P1_vec++;
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_spread2(zero, in, z);
            GFNI_MACC2(t, z, B1 + (c / 2) * KP);
            c++;
        }
        for (; c + 1 < V_MAX; c += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P1, P1_vec), a);
            gfni_load_mvec(PINV(P1, P1_vec + 1), b);
            P1_vec += 2;
            gfni_spread2(a, b, z);
            GFNI_MACC2(t, z, B1 + (c / 2) * KP);
        }
        if (c < V_MAX) {        /* odd V_MAX: a last column with no partner */
            __m512i in[GFNI_R], z[ZP];
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_load_mvec(PINV(P1, P1_vec), in);
            P1_vec++;
            gfni_spread2(in, zero, z);
            GFNI_MACC2(t, z, B1 + (c / 2) * KP);
        }

        // P2 * S2
        size_t cc = 0;
        for (; cc + 1 < O_MAX; cc += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(PINV(P2, r * O_MAX + cc), a);
            gfni_load_mvec(PINV(P2, r * O_MAX + cc + 1), b);
            gfni_spread2(a, b, z);
            GFNI_MACC2(t, z, B2 + (cc / 2) * KP);
        }
        if (cc < O_MAX) {       /* odd O_MAX */
            __m512i in[GFNI_R], z[ZP];
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_load_mvec(PINV(P2, r * O_MAX + cc), in);
            gfni_spread2(in, zero, z);
            GFNI_MACC2(t, z, B2 + (cc / 2) * KP);
        }

        for (size_t lp = 0; lp < KP; lp++) {
            uint64_t *q = PS2X2_PAIR(acc, r, lp);
            for (size_t w = 0; w < ZP; w++) _mm512_storeu_si512((void *) (q + 8 * w), t[lp][w]);
        }
    }
}

// P3*S2 -> P3: o x o, S2: o x k // P3 upper triangular
static
inline void P3_times_S2(const uint64_t *P3, const uint64_t *S2_tabs, uint64_t *acc) {
    uint64_t B[((O_MAX + 1) / 2) * KP];
    GFNI_BLK_TABLE(B, S2_tabs, S2_TAB_STRIDE, O_MAX);
    size_t cols_used = 0;
    for (size_t r = 0; r < O_MAX; r++) {
        __m512i t[KP][ZP];
        for (size_t lp = 0; lp < KP; lp++) {
            for (size_t q = 0; q < ZP; q++) t[lp][q] = _mm512_setzero_si512();
        }
        size_t c = r;
        if (c & 1) {
            __m512i in[GFNI_R], z[ZP];
            gfni_load_mvec(P3 + cols_used, in);
            cols_used += M_VEC_LIMBS_MAX;
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_spread2(zero, in, z);
            GFNI_MACC2(t, z, B + (c / 2) * KP);
            c++;
        }
        for (; c + 1 < O_MAX; c += 2) {
            __m512i a[GFNI_R], b[GFNI_R], z[ZP];
            gfni_load_mvec(P3 + cols_used, a);
            cols_used += M_VEC_LIMBS_MAX;
            gfni_load_mvec(P3 + cols_used, b);
            cols_used += M_VEC_LIMBS_MAX;
            gfni_spread2(a, b, z);
            GFNI_MACC2(t, z, B + (c / 2) * KP);
        }
        if (c < O_MAX) {        /* odd O_MAX */
            __m512i in[GFNI_R], z[ZP];
            __m512i zero[GFNI_R];
            for (size_t q = 0; q < GFNI_R; q++) zero[q] = _mm512_setzero_si512();
            gfni_load_mvec(P3 + cols_used, in);
            cols_used += M_VEC_LIMBS_MAX;
            gfni_spread2(in, zero, z);
            GFNI_MACC2(t, z, B + (c / 2) * KP);
        }
        for (size_t lp = 0; lp < KP; lp++) {
            uint64_t *q = PS2X2_PAIR(acc, r, lp);
            for (size_t w = 0; w < ZP; w++) _mm512_storeu_si512((void *) (q + 8 * w), t[lp][w]);
        }
    }
}

static inline
void compute_M_and_VPV(const mayo_params_t *p, const unsigned char *Vdec, const uint64_t *L, const uint64_t *P1,
                       uint64_t *VL, uint64_t *VP1V) {
    (void) p;
    alignas(MAYO_MULTAB_ALIGN) uint64_t V_tabs[MAYO_V_MULTABS_N];
    mayo_V_multabs(Vdec, V_tabs);

    // M
    Vt_times_L(L, V_tabs, VL);

    // VP1V
    uint64_t Pv[V_MAX * K_MAX * M_VEC_LIMBS_MAX + MAYO_MVEC_SLACK] = {0};
    P1_times_Vt(P1, V_tabs, Pv);
    Vt_times_Pv(Pv, V_tabs, VP1V);

    mayo_secure_clear(V_tabs, sizeof(V_tabs));
    mayo_secure_clear(Pv, sizeof(Pv));
}

static inline
void compute_P3(const mayo_params_t *p, const uint64_t *P1, uint64_t *P2, const unsigned char *O, uint64_t *P3) {
    (void) p;
    uint64_t O_tabs[MAYO_O_MULTABS_N];
    mayo_O_multabs(O, O_tabs);
    P1_times_O(P1, O_tabs, P2);
    Ot_times_P1O_P2(P2, O_tabs, P3);
    mayo_secure_clear(O_tabs, sizeof(O_tabs));
}

// compute P * S^t = [ P1  P2 ] * [S1] = [P1*S1 + P2*S2]
//                   [  0  P3 ]   [S2]   [        P3*S2]
// compute S * PS  = [ S1 S2 ] * [ P1*S1 + P2*S2 = P1 ] = [ S1*P1 + S2*P2 ]
//                               [         P3*S2 = P2 ]
static inline void m_calculate_PS_SPS(const mayo_params_t *p, const uint64_t *P1, const uint64_t *P2, const uint64_t *P3,
                                      const unsigned char *S, const uint64_t *lambda, uint64_t *SPS) {
    (void) p;
    const int o = PARAM_NAME(o);
    const int v = PARAM_NAME(v);
    const int k = PARAM_NAME(k);
    const int n = o + v;
    unsigned char S1[V_MAX * K_MAX]; // == N-O, K
    unsigned char S2[O_MAX * K_MAX]; // == O, K
    unsigned char *s1_write = S1;
    unsigned char *s2_write = S2;

    for (int r = 0; r < k; r++) {
        for (int c = 0; c < n; c++) {
            if (c < v) {
                *(s1_write++) = S[r * n + c];
            } else {
                *(s2_write++) = S[r * n + c];
            }
        }
    }

    alignas(64) uint64_t PS[N_MAX * KP * ZP * 8]; // interleaved pairs, fully written below
#define MAYO_PS_ROWS_OFF (V_MAX * KP * ZP * 8)
    MAYO_ZERO_SLACK(PS, (size_t) N_MAX * K_MAX * M_VEC_LIMBS_MAX);

    alignas(MAYO_MULTAB_ALIGN) uint64_t S1_tabs[MAYO_V_MULTABS_N];
    alignas(MAYO_MULTAB_ALIGN) uint64_t S2_tabs[MAYO_S2_MULTABS_N];
    mayo_S1_multabs(S1, S1_tabs);
    mayo_S2_multabs(S2, S2_tabs);

    P1_times_S1_plus_P2_times_S2(P1, P2, S1_tabs, S2_tabs, PS);
    P3_times_S2(P3, S2_tabs, PS + MAYO_PS_ROWS_OFF); // upper triangular

    /* lambda occupies column 0 of PS, which in this layout is the even half of column
       pair 0; spread it into that half and xor it into the pair slot, the same full-width
       slot the producers above store. The S^t * PS pass then applies S[i][r] to it, so
       SPS[i][0] picks up Lambda(s_i) for n m-vector XORs and no multiplications. */
    {
        __m512i zero[GFNI_R];
        for (size_t q = 0; q < GFNI_R; q++) {
            zero[q] = _mm512_setzero_si512();
        }
        for (int r = 0; r < n; r++) {
            uint64_t *slot = (r < v) ? PS2X2_PAIR(PS, r, 0)
                                     : PS2X2_PAIR(PS + MAYO_PS_ROWS_OFF, r - v, 0);
            __m512i in[GFNI_R], z[ZP];
            gfni_load_mvec(lambda + (size_t) r * M_VEC_LIMBS_MAX, in);
            gfni_spread2(in, zero, z);
            for (size_t w = 0; w < ZP; w++) {
                _mm512_storeu_si512((void *) (slot + 8 * w),
                                    _mm512_xor_si512(_mm512_loadu_si512((const void *) (slot + 8 * w)), z[w]));
            }
        }
    }

    // S^T * PS = S1^t*PS1 + S2^t*PS2
    S1t_times_PS1(PS, S1_tabs, SPS);
    S2t_times_PS2(PS + MAYO_PS_ROWS_OFF, S2_tabs, SPS);
}


/* 8x8 transpose of 64-bit lanes: dst row i = column i of the input rows */
static inline void gfni_transpose_8x8_epi64(__m512i *r) {
    __m512i t[8], u[8];
    t[0] = _mm512_unpacklo_epi64(r[0], r[1]);
    t[1] = _mm512_unpackhi_epi64(r[0], r[1]);
    t[2] = _mm512_unpacklo_epi64(r[2], r[3]);
    t[3] = _mm512_unpackhi_epi64(r[2], r[3]);
    t[4] = _mm512_unpacklo_epi64(r[4], r[5]);
    t[5] = _mm512_unpackhi_epi64(r[4], r[5]);
    t[6] = _mm512_unpacklo_epi64(r[6], r[7]);
    t[7] = _mm512_unpackhi_epi64(r[6], r[7]);
    u[0] = _mm512_shuffle_i64x2(t[0], t[2], 0x88);
    u[1] = _mm512_shuffle_i64x2(t[4], t[6], 0x88);
    u[2] = _mm512_shuffle_i64x2(t[1], t[3], 0x88);
    u[3] = _mm512_shuffle_i64x2(t[5], t[7], 0x88);
    u[4] = _mm512_shuffle_i64x2(t[0], t[2], 0xdd);
    u[5] = _mm512_shuffle_i64x2(t[4], t[6], 0xdd);
    u[6] = _mm512_shuffle_i64x2(t[1], t[3], 0xdd);
    u[7] = _mm512_shuffle_i64x2(t[5], t[7], 0xdd);
    r[0] = _mm512_shuffle_i64x2(u[0], u[1], 0x88);
    r[4] = _mm512_shuffle_i64x2(u[0], u[1], 0xdd);
    r[1] = _mm512_shuffle_i64x2(u[2], u[3], 0x88);
    r[5] = _mm512_shuffle_i64x2(u[2], u[3], 0xdd);
    r[2] = _mm512_shuffle_i64x2(u[4], u[5], 0x88);
    r[6] = _mm512_shuffle_i64x2(u[4], u[5], 0xdd);
    r[3] = _mm512_shuffle_i64x2(u[6], u[7], 0x88);
    r[7] = _mm512_shuffle_i64x2(u[6], u[7], 0xdd);
}

static inline void gfni_transpose_4x4_epi64(__m256i *r) {
    __m256i t0 = _mm256_unpacklo_epi64(r[0], r[1]);
    __m256i t1 = _mm256_unpackhi_epi64(r[0], r[1]);
    __m256i t2 = _mm256_unpacklo_epi64(r[2], r[3]);
    __m256i t3 = _mm256_unpackhi_epi64(r[2], r[3]);
    r[0] = _mm256_permute2x128_si256(t0, t2, 0x20);
    r[1] = _mm256_permute2x128_si256(t1, t3, 0x20);
    r[2] = _mm256_permute2x128_si256(t0, t2, 0x31);
    r[3] = _mm256_permute2x128_si256(t1, t3, 0x31);
}

/* see the fallback in arithmetic.h. Eight columns are loaded whole and transposed in
   registers rather than gathered. The loads run past the last column into Mtmp's
   MAYO_MVEC_SLACK tail; every index and the shift amount are public. */
/* dst[i] ^= each nibble of src[i] times coeff, one affine multiply per eight words */
#define MAYO_HAVE_GF16_MULC_XOR 1
static inline void m_gf16_mulc_xor(uint64_t *dst, const uint64_t *src, int n, unsigned char coeff, const unsigned char *t16) {
    (void) t16;   /* the affine matrix folds away for a constant coeff */
    const __m512i mat = _mm512_set1_epi64((long long) mayo_gfni_tab(coeff));
    const int batched = n & ~7;
    int i = 0;
    for (; i < batched; i += 8) {
        __m512i v = _mm512_loadu_si512((const void *) (src + i));
        _mm512_storeu_si512((void *) (dst + i),
            _mm512_xor_si512(_mm512_loadu_si512((const void *) (dst + i)),
                             _mm512_gf2p8affine_epi64_epi8(v, mat, 0)));
    }
    if (i < n) {
        const __mmask8 km = (__mmask8) ((1u << (n - i)) - 1u);
        __m512i v = _mm512_maskz_loadu_epi64(km, src + i);
        __m512i d = _mm512_maskz_loadu_epi64(km, dst + i);
        _mm512_mask_storeu_epi64(dst + i, km,
            _mm512_xor_si512(d, _mm512_gf2p8affine_epi64_epi8(v, mat, 0)));
    }
}

#define MAYO_HAVE_COMPUTE_A_ACC 1
static inline void m_compute_A_acc(uint64_t *A, size_t A_width, int row0, int col,
                                   const uint64_t *M, int o, int mvl, int bits) {
    const __m128i shl = _mm_cvtsi32_si128(bits);
    const __m128i shr = _mm_cvtsi32_si128(64 - bits);
    int c = 0;
    /* whole 8-column blocks only: a masked store costs a full store either way, so a
       partial block moves a fraction of the data for the same price */
    for (; c + 8 <= o; c += 8) {
        for (int l0 = 0; l0 < mvl; l0 += 8) {
            const int lv = (mvl - l0) < 8 ? (mvl - l0) : 8;
            __m512i r[8];
            for (int t = 0; t < 8; t++) {
                r[t] = _mm512_loadu_si512((const void *) (M + (size_t) (c + t) * mvl + l0));
            }
            gfni_transpose_8x8_epi64(r);
            for (int l = 0; l < lv; l++) {
                uint64_t *d = A + (size_t) (row0 + l0 + l) * A_width + col + c;
                _mm512_storeu_si512((void *) d,
                    _mm512_xor_si512(_mm512_loadu_si512((const void *) d), _mm512_sll_epi64(r[l], shl)));
                if (bits > 0) {
                    uint64_t *d1 = A + (size_t) (row0 + l0 + l + 1) * A_width + col + c;
                    _mm512_storeu_si512((void *) d1,
                        _mm512_xor_si512(_mm512_loadu_si512((const void *) d1), _mm512_srl_epi64(r[l], shr)));
                }
            }
        }
    }
    /* the remainder in narrower whole blocks, same idea: 4 then 2 columns, all with
       plain stores. MAYO-3 leaves exactly 2 and MAYO-5 exactly 4. */
    for (; c + 4 <= o; c += 4) {
        for (int l0 = 0; l0 < mvl; l0 += 4) {
            const int lv = (mvl - l0) < 4 ? (mvl - l0) : 4;
            __m256i r[4];
            for (int t = 0; t < 4; t++) {
                r[t] = _mm256_loadu_si256((const __m256i *) (M + (size_t) (c + t) * mvl + l0));
            }
            gfni_transpose_4x4_epi64(r);
            for (int l = 0; l < lv; l++) {
                uint64_t *d = A + (size_t) (row0 + l0 + l) * A_width + col + c;
                _mm256_storeu_si256((__m256i *) d,
                    _mm256_xor_si256(_mm256_loadu_si256((const __m256i *) d), _mm256_sll_epi64(r[l], shl)));
                if (bits > 0) {
                    uint64_t *d1 = A + (size_t) (row0 + l0 + l + 1) * A_width + col + c;
                    _mm256_storeu_si256((__m256i *) d1,
                        _mm256_xor_si256(_mm256_loadu_si256((const __m256i *) d1), _mm256_srl_epi64(r[l], shr)));
                }
            }
        }
    }
    for (; c + 2 <= o; c += 2) {
        for (int l0 = 0; l0 < mvl; l0 += 2) {
            const int lv = (mvl - l0) < 2 ? (mvl - l0) : 2;
            __m128i a0 = _mm_loadu_si128((const __m128i *) (M + (size_t) c * mvl + l0));
            __m128i a1 = _mm_loadu_si128((const __m128i *) (M + (size_t) (c + 1) * mvl + l0));
            __m128i r0 = _mm_unpacklo_epi64(a0, a1);
            __m128i r1 = _mm_unpackhi_epi64(a0, a1);
            const __m128i rr[2] = { r0, r1 };
            for (int l = 0; l < lv; l++) {
                uint64_t *d = A + (size_t) (row0 + l0 + l) * A_width + col + c;
                _mm_storeu_si128((__m128i *) d,
                    _mm_xor_si128(_mm_loadu_si128((const __m128i *) d), _mm_sll_epi64(rr[l], shl)));
                if (bits > 0) {
                    uint64_t *d1 = A + (size_t) (row0 + l0 + l + 1) * A_width + col + c;
                    _mm_storeu_si128((__m128i *) d1,
                        _mm_xor_si128(_mm_loadu_si128((const __m128i *) d1), _mm_srl_epi64(rr[l], shr)));
                }
            }
        }
    }
    if (c < o) {    /* the narrowest block leaves at most one column */
        for (int l = 0; l < mvl; l++) {
            const uint64_t v = M[(size_t) c * mvl + l];
            A[(size_t) (row0 + l) * A_width + col + c] ^= v << bits;
            if (bits > 0) {
                A[(size_t) (row0 + l + 1) * A_width + col + c] ^= v >> (64 - bits);
            }
        }
    }
}

#undef K_OVER_2
#endif
