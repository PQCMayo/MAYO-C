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
#include <m_vec_arithmetic.h>

// GF(16) arithmetic using GFNI, on zmm with AVX-512, else on ymm. Multiplication by a constant is an
// F2-linear map on the four nibble bits, i.e. an 8x8 bit-matrix that
// vgf2p8affineqb applies to a whole register of packed bytes at once, so m-vecs stay packed
// (see the width layer below) and a multiply-accumulate is one affine
// plus one xor.

// the first producer into VKtmp/SPS assigns, so mayo.c need not pre-zero them
#define MAYO_ACC_STORES_FIRST 1
// P1/P2 are read straight from the packed keystream layout (m/2 bytes per m-vec),
// skipping the unpack spread. PINV turns an m-vec index into the packed byte
// address -- the identity when m % 16 == 0. A load reads M_VEC_LIMBS_MAX limbs,
// so its tail nibbles are next-vec pad (harmless).
#define MAYO_PACKED_P1P2 1
#if (M_MAX % 16) == 0
#define PINV(base, vec)      ((const uint64_t *)(base) + (size_t)(vec) * M_VEC_LIMBS_MAX)
#else
#define PINV(base, vec)      ((const uint64_t *)((const uint8_t *)(base) + (size_t)(vec) * (size_t)(M_MAX / 2)))
#endif

typedef uint64_t mayo_multab_t;
/* V's "multabs" are its 2x2 block table, built once per attempt for all four V kernels */
#define MAYO_V_MULTABS_N  (((V_MAX + 1) / 2) * ((K_MAX + 1) / 2))
#define MAYO_S1_MULTABS_N (K_MAX * V_MAX)
#define MAYO_S2_MULTABS_N (K_MAX * O_MAX)
#define MAYO_O_MULTABS_N (O_MAX * V_MAX)
// select the multab sign orchestration in arithmetic.h
#define MAYO_HAVE_MULTAB_SIGN 1

// arithmetic.h's sign orchestration builds its tables through these names
#define mayo_O_multabs mayo_gfni_O_multabs
#define mayo_V_multabs mayo_gfni_V_multabs

// The per-element tables (mayo_gfni_tabs) are a straight elementwise map, so V/S1/S2's
// stay in their input (row-major) order and the block builders walk them with a stride
// of V_TAB_STRIDE/S2_TAB_STRIDE.
#define V_TAB_STRIDE  V_MAX
#define S2_TAB_STRIDE O_MAX

static
inline void mayo_gfni_O_multabs(const unsigned char *O, uint64_t *O_tabs) {
    mayo_gfni_tabs(O, O_tabs, (size_t) V_MAX * O_MAX);
}


// Vector width. The kernels are written once against gfv: zmm when AVX-512 is there,
// otherwise ymm (GFNI's VEX form, e.g. Alder/Meteor Lake, Zen 4/5 at 256 bits).
#ifdef MAYO_GFNI_ZMM
typedef __m512i gfv;
#define GFNI_W 8                /* limbs per register */
#define GFNI_ACCS 24            /* accumulator budget, 32 registers */
#define gfv_zero() _mm512_setzero_si512()
#define gfv_set1_64(x) _mm512_set1_epi64((long long) (x))
#define gfv_set1_8(x) _mm512_set1_epi8(x)
#define gfv_affine(x, m) _mm512_gf2p8affine_epi64_epi8((x), (m), 0)
#define gfv_loadu(p) _mm512_loadu_si512((const void *) (p))
#define gfv_storeu(p, x) _mm512_storeu_si512((void *) (p), (x))
#define gfv_xor(a, b) _mm512_xor_si512((a), (b))
#define gfv_xor3(a, b, c) _mm512_ternarylogic_epi64((a), (b), (c), 0x96)   /* one uop, not two */
/* nibble moves for the (de)spread, each byte on its own; the other nibble of the result
   is unspecified, every use masks it off. nib_down: high nibble to low; nib_up: low to
   high. With VBMI a vpmultishiftqb (shuffle port) instead of a 512-bit shift (port 0,
   which the affines need): per qword, ctl byte i picks the 8 bits from bit ctl[i]. */
#ifdef __AVX512VBMI__
#define GFNI_MS_ID 0x3830282018100800LL
#define GFNI_MS_DOWN 0x3C342C241C140C04LL
#define GFNI_MS_UP 0x342C241C140C043CLL
#define gfv_nib_down(x) _mm512_multishift_epi64_epi8(_mm512_set1_epi64(GFNI_MS_DOWN), (x))
#define gfv_nib_up(x) _mm512_multishift_epi64_epi8(_mm512_set1_epi64(GFNI_MS_UP), (x))
/* the same in the high (resp. low) 256 bits only, the identity in the other */
#define gfv_nib_down_hi(x) _mm512_multishift_epi64_epi8(_mm512_setr_epi64(GFNI_MS_ID, GFNI_MS_ID, \
        GFNI_MS_ID, GFNI_MS_ID, GFNI_MS_DOWN, GFNI_MS_DOWN, GFNI_MS_DOWN, GFNI_MS_DOWN), (x))
#define gfv_nib_up_lo(x) _mm512_multishift_epi64_epi8(_mm512_setr_epi64(GFNI_MS_UP, GFNI_MS_UP, \
        GFNI_MS_UP, GFNI_MS_UP, GFNI_MS_ID, GFNI_MS_ID, GFNI_MS_ID, GFNI_MS_ID), (x))
#else
#define gfv_nib_down(x) _mm512_srli_epi16((x), 4)
#define gfv_nib_up(x) _mm512_slli_epi16((x), 4)
#define gfv_nib_down_hi(x) _mm512_mask_srli_epi16((x), 0xFFFF0000u, (x), 4)
#define gfv_nib_up_lo(x) _mm512_mask_slli_epi16((x), 0x0000FFFFu, (x), 4)
#endif
/* (lo & a) | (~lo & b) */
#define gfv_sel(lo, a, b) _mm512_ternarylogic_epi64((lo), (a), (b), 0xCA)
/* cat_lo: the low halves of a and b; dup_hi: the high half of a, twice */
#define gfv_cat_lo(a, b) _mm512_shuffle_i32x4((a), (b), 0x44)
#define gfv_dup_hi(a) _mm512_shuffle_i32x4((a), (a), 0xEE)
#else
typedef __m256i gfv;
#define GFNI_W 4
#define GFNI_ACCS 12            /* 16 registers: 12 accumulators, the pair and a matrix */
#define gfv_zero() _mm256_setzero_si256()
#define gfv_set1_64(x) _mm256_set1_epi64x((long long) (x))
#define gfv_set1_8(x) _mm256_set1_epi8(x)
#define gfv_affine(x, m) _mm256_gf2p8affine_epi64_epi8((x), (m), 0)
#define gfv_loadu(p) _mm256_loadu_si256((const __m256i *) (p))
#define gfv_storeu(p, x) _mm256_storeu_si256((__m256i *) (p), (x))
#define gfv_xor(a, b) _mm256_xor_si256((a), (b))
#define gfv_xor3(a, b, c) gfv_xor(gfv_xor((a), (b)), (c))
#define gfv_nib_down(x) _mm256_srli_epi16((x), 4)
#define gfv_nib_up(x) _mm256_slli_epi16((x), 4)
#define gfv_sel(lo, a, b) _mm256_or_si256(_mm256_and_si256((lo), (a)), _mm256_andnot_si256((lo), (b)))
#define gfv_cat_lo(a, b) _mm256_permute2x128_si256((a), (b), 0x20)
#define gfv_dup_hi(a) _mm256_permute2x128_si256((a), (a), 0x11)
#endif

/* The m-vec is swept in pieces, each kernel running once per piece, so that a piece's
   accumulators fit the register file: zmm takes the whole m-vec in one piece (two
   registers for m > 128), ymm one register per piece. (Half-register ymm pieces, one
   accumulator per output pair, measured 5-15% slower: twice the matrix broadcasts and
   spreads per affine outweigh the spills they avoid.) The piece loop is unrolled, so
   pc is a constant in each copy of the body; with one piece it folds away entirely. */
#ifdef MAYO_GFNI_ZMM
#define GFNI_RP ((M_VEC_LIMBS_MAX + GFNI_W - 1) / GFNI_W)   /* registers per piece */
#else
#define GFNI_RP 1
#endif
#define GFNI_PCL (GFNI_RP * GFNI_W)                          /* limbs per piece */
#define GFNI_NPC ((M_VEC_LIMBS_MAX + GFNI_PCL - 1) / GFNI_PCL)
/* at most one per function. One piece is no loop at all: gcc schedules the kernels
   differently inside even a one-trip loop. */
#define GFNI_PRAGMA_(x) _Pragma(#x)
#if GFNI_NPC == 1
#define GFNI_FOR_PIECES(pc) const size_t pc = 0; (void) pc; if (1)
#elif defined(__clang__)
#define GFNI_FOR_PIECES(pc) GFNI_PRAGMA_(clang loop unroll(full)) for (size_t pc = 0; pc < GFNI_NPC; pc++)
#else
#define GFNI_FOR_PIECES(pc) GFNI_PRAGMA_(GCC unroll 4) for (size_t pc = 0; pc < GFNI_NPC; pc++)
#endif

// Loads are full-width even though an m-vec is only M_VEC_LIMBS_MAX limbs, and the
// stores are masked to its limbs, so the extra lanes are discarded (gcc schedules masked
// loads badly: 2x slower on the products). Buffers carry MAYO_MVEC_SLACK limbs so the
// wider read stays in bounds.
#ifdef MAYO_GFNI_ZMM
/* one piece: the m-vec is one zmm, or two for m > 128, the last one partial */
#define GFNI_TAIL_OFF (8 * (GFNI_RP - 1))
#define GFNI_TAIL_MASK ((__mmask8)((1u << (M_VEC_LIMBS_MAX - GFNI_TAIL_OFF)) - 1))
/* a last register of 2..4 limbs is loaded into both 256-bit halves, so gfni_spread2 builds
   its merged pair register in place instead of shuffling two halves together */
#define GFNI_HALVES (M_VEC_LIMBS_MAX - GFNI_TAIL_OFF > 1 && M_VEC_LIMBS_MAX - GFNI_TAIL_OFF <= 4)
static inline void gfni_load_mvec(size_t pc, const uint64_t *p, gfv in[GFNI_RP]) {
    (void) pc;
#if GFNI_RP == 2
    in[0] = _mm512_loadu_si512((const void *) p);
#endif
#if M_VEC_LIMBS_MAX % 8 == 1 && M_VEC_LIMBS_MAX > 1
    /* a 1-limb last register: just the limb, so its spread leaves qwords 1 and 5 zero for
       the quad steps (see GFNI_LAST1) */
    in[GFNI_RP - 1] = _mm512_zextsi128_si512(_mm_loadl_epi64((const __m128i *) (p + GFNI_TAIL_OFF)));
#elif GFNI_HALVES
    in[GFNI_RP - 1] = _mm512_broadcast_i64x4(_mm256_loadu_si256((const __m256i *) (p + GFNI_TAIL_OFF)));
#else
    in[GFNI_RP - 1] = _mm512_loadu_si512((const void *) (p + GFNI_TAIL_OFF));
#endif
}
/* at m = 64 a despread pair fills one register exactly: one store, not two masked halves */
#if M_VEC_LIMBS_MAX == 4
#define GFNI_PAIR_STORE_FITS 1
static inline void gfni_xor_store_pair(size_t pc, uint64_t *p, const gfv x[GFNI_RP], const gfv y[GFNI_RP]) {
    (void) pc;
    const gfv v = _mm512_shuffle_i32x4(x[0], y[0], 0x44);
    _mm512_storeu_si512((void *) p, _mm512_xor_si512(_mm512_loadu_si512((const void *) p), v));
}
#endif
// acc ^= t, touching exactly M_VEC_LIMBS_MAX limbs
static inline void gfni_xor_store_mvec(size_t pc, uint64_t *p, const gfv t[GFNI_RP]) {
    (void) pc;
#if GFNI_RP == 2
    _mm512_storeu_si512((void *) p, _mm512_xor_si512(_mm512_loadu_si512((const void *) p), t[0]));
#endif
    _mm512_mask_storeu_epi64((void *) (p + GFNI_TAIL_OFF), GFNI_TAIL_MASK,
                             _mm512_xor_si512(_mm512_loadu_si512((const void *) (p + GFNI_TAIL_OFF)),
                                              t[GFNI_RP - 1]));
}
// acc = t (no accumulate). Producers that write every element of their output exactly
// once use this, so the target needs no zeroing beforehand; it also drops a load.
static inline void gfni_store_mvec(size_t pc, uint64_t *p, const gfv t[GFNI_RP]) {
    (void) pc;
#if GFNI_RP == 2
    _mm512_storeu_si512((void *) p, t[0]);
#endif
    _mm512_mask_storeu_epi64((void *) (p + GFNI_TAIL_OFF), GFNI_TAIL_MASK, t[GFNI_RP - 1]);
}
#else
/* GFNI_PCL limbs of the m-vec, fewer in the last piece. A load reads a whole ymm; the
   lanes past the piece are dropped by the spread or the store. */
#define GFNI_PC_LIMBS(pc) (M_VEC_LIMBS_MAX - GFNI_PCL * (pc) < GFNI_PCL ? M_VEC_LIMBS_MAX - GFNI_PCL * (pc) : GFNI_PCL)
static inline void gfni_load_mvec(size_t pc, const uint64_t *p, gfv in[GFNI_RP]) {
    if (GFNI_PC_LIMBS(pc) == 1) {   /* alone: its spread pair leaves qwords 1 and 3 zero (GFNI_QUAD) */
        in[0] = _mm256_zextsi128_si256(_mm_loadl_epi64((const __m128i *) (p + GFNI_PCL * pc)));
    } else {
        in[0] = _mm256_loadu_si256((const __m256i *) (p + GFNI_PCL * pc));
    }
}
/* p[0 .. n) = x[0 .. n); n is a constant per piece */
static inline void gfni_store_part(uint64_t *p, gfv x, size_t n) {
    if (n == 4) {
        _mm256_storeu_si256((__m256i *) p, x);
    } else if (n == 1) {
        _mm_storel_epi64((__m128i *) p, _mm256_castsi256_si128(x));
    } else {
        const __m256i m = _mm256_cmpgt_epi64(_mm256_set1_epi64x((long long) n), _mm256_setr_epi64x(0, 1, 2, 3));
        _mm256_maskstore_epi64((long long *) p, m, x);
    }
}
// acc ^= t, touching exactly the piece's limbs
static inline void gfni_xor_store_mvec(size_t pc, uint64_t *p, const gfv t[GFNI_RP]) {
    p += GFNI_PCL * pc;
    gfni_store_part(p, _mm256_loadu_si256((const __m256i *) p) ^ t[0], GFNI_PC_LIMBS(pc));
}
// acc = t (no accumulate), as above
static inline void gfni_store_mvec(size_t pc, uint64_t *p, const gfv t[GFNI_RP]) {
    gfni_store_part(p + GFNI_PCL * pc, t[0], GFNI_PC_LIMBS(pc));
}
#endif

// 2x2 GF(16) blocks: one affine can apply any F2-linear map to a byte, so with two
// m-vecs nibble-interleaved (byte j = x's nibble j low, y's high) it computes
// (a*x + b*y, c*x + d*y) -- two inputs and two outputs per affine instead of one.
// Byte j of the matrix qword is output bit 7-j, so out-low sits in bytes 4..7.

// slots, not pairs: an odd k leaves a last slot whose second row is zero, so its
// block computes the leftover output in the low half and nothing in the high
#define KP ((K_MAX + 1) / 2)
/* registers an interleaved pair spans: one byte per GF(16) element of the two m-vecs.
   A ymm piece's pair takes two registers, one if its limbs fit half a register. */
#ifdef MAYO_GFNI_ZMM
#define GFNI_ZPP(pc) ((M_MAX + 63) / 64)
#else
#define GFNI_ZPP(pc) (GFNI_PC_LIMBS(pc) > 2 ? 2 : 1)
#endif
/* the first piece is never the short one, so it spans the most */
#define GFNI_ZPMAX GFNI_ZPP(0)
/* a stored pair: one register per half register of m-vec, the pieces back to back */
#define GFNI_SLOT (GFNI_W * ((M_VEC_LIMBS_MAX + GFNI_W / 2 - 1) / (GFNI_W / 2)))
#define GFNI_ZOFF(pc, w) (GFNI_W * (GFNI_ZPMAX * (pc) + (w)))
// the consumers hold one accumulator per output and register, so chunk the outputs to
// keep that inside the register file
#define OCH ((GFNI_ACCS / GFNI_ZPMAX) < O_MAX ? (GFNI_ACCS / GFNI_ZPMAX) : O_MAX)
#define GFNI_2X2_AM 0x0F0F0F0F00000000ULL
#define GFNI_2X2_BM 0x00000000F0F0F0F0ULL
static inline uint64_t gfni_blk2(uint64_t ta, uint64_t tb, uint64_t tc, uint64_t td) {
    return (ta & GFNI_2X2_AM) | ((tb & GFNI_2X2_AM) << 4)
         | ((tc & GFNI_2X2_BM) >> 4) | (td & GFNI_2X2_BM);
}
/* Quad steps: four source m-vecs (two spread pairs) per step, against the blocks BA
   (pair a, b) and BB (pair c, d) of NP output pairs. */
#ifdef MAYO_GFNI_ZMM
/* zmm: every set takes them, each accumulator updated once for both pairs (one three-way
   xor). A 1-limb last register (MAYO_5, GFNI_LAST1) is tail-packed, see below. */
#define GFNI_QUAD(pc) 1
#define GFNI_LAST1 (M_VEC_LIMBS_MAX % GFNI_W == 1 && M_VEC_LIMBS_MAX > 1)
#define GFNI_TAILPACK GFNI_LAST1
#define GFNI_FOLD(z) (z)   /* only the ymm quad steps leave a split sum */
/* Tail packing: the last register of a spread pair is then one qword of even and one of
   odd elements, so it is spread out to [ab, ab, cd, cd | ab', ab', cd', cd'] and output
   pairs (p, p+1) share one accumulator t[p][last] (p even) with the four matrices
   [M_ab(p), M_ab(p+1), M_cd(p), M_cd(p+1)] per half: two input pairs by two output pairs
   per affine. GFNI_TAIL_UNPACK restores the plain layout before the despread. */
/* the m-vec's last limb; packed P1 m-vecs need not be 8-byte aligned (PINV) */
static inline long long gfni_tail_limb(const void *p) {
    long long v;
    memcpy(&v, (const unsigned char *) p + 8 * GFNI_TAIL_OFF, sizeof(v));
    return v;
}
static inline gfv gfni_spread4_last(const uint64_t *a, const uint64_t *b, const uint64_t *c, const uint64_t *d) {
    const long long la = gfni_tail_limb(a), lb = gfni_tail_limb(b), lc = gfni_tail_limb(c), ld = gfni_tail_limb(d);
    const __m512i x = _mm512_setr_epi64(la, la, lc, lc, la, la, lc, lc);
    const __m512i y = _mm512_setr_epi64(lb, lb, ld, ld, lb, lb, ld, ld);
    return gfv_sel(gfv_set1_8(0x0f), _mm512_srlv_epi64(x, _mm512_setr_epi64(0, 0, 0, 0, 4, 4, 4, 4)),
                   _mm512_sllv_epi64(y, _mm512_setr_epi64(4, 4, 4, 4, 0, 0, 0, 0)));
}
#define GFNI_QUAD_STEP(t, a, b, c, d, BA, BB, NP)                                         \
    do {                                                                                     \
        gfv _xa[GFNI_RP], _xb[GFNI_RP], _xc[GFNI_RP], _xd[GFNI_RP], _za[GFNI_ZPMAX], _zc[GFNI_ZPMAX]; \
        gfni_load_mvec(pc, (a), _xa); gfni_load_mvec(pc, (b), _xb);                          \
        gfni_load_mvec(pc, (c), _xc); gfni_load_mvec(pc, (d), _xd);                          \
        gfni_spread2(pc, _xa, _xb, _za);                                                     \
        gfni_spread2(pc, _xc, _xd, _zc);                                                     \
        const gfv _zl = GFNI_LAST1 ? gfni_spread4_last((a), (b), (c), (d)) : gfv_zero();     \
        GFNI_PRAGMA_(GCC unroll 16)                                                          \
        for (size_t _p = 0; _p < (size_t)(NP); _p++) {                                       \
            gfv _ma = gfv_set1_64((BA)[_p]), _mc = gfv_set1_64((BB)[_p]);                                       \
            __asm__("" : "+v"(_ma), "+v"(_mc));                                              \
            for (size_t _q = 0; _q + GFNI_LAST1 < GFNI_ZPP(pc); _q++) {                       \
                (t)[_p][_q] = gfv_xor3((t)[_p][_q], gfv_affine(_za[_q], _ma), gfv_affine(_zc[_q], _mc)); \
            }                                                                                \
        }                                                                                    \
        for (size_t _p = 0; GFNI_LAST1 && _p < (size_t)(NP); _p += 2) {                      \
            const gfv _m = GFNI_TAIL_MAT4((BA), (BB), _p);                               \
            (t)[_p][GFNI_ZPP(pc) - 1] ^= gfv_affine(_zl, _m);                                \
        }                                                                                    \
    } while (0)
/* [M_a(p), M_a(p+1), M_c(p), M_c(p+1)] in both halves: two 16-byte loads from the block
   tables (adjacent output pairs are adjacent entries; MAYO_5 has an even KP and OP) */
#define GFNI_TAIL_MAT4(BA, BB, p)                                                        \
    _mm512_broadcast_i64x4(_mm256_loadu2_m128i((const __m128i *) ((BB) + (p)),               \
                                               (const __m128i *) ((BA) + (p))))
/* a pair step's tail: z's last register is [ab, 0, 0, 0 | ab', 0, 0, 0] (1-limb load) */
#define GFNI_TAIL_PAIR(t, zl, BLK, NP)                                                       \
    do {                                                                                     \
        const gfv _zl2 = _mm512_permutex_epi64((zl), 0x50);   /* [ab, ab, 0, 0 | ...] */      \
        for (size_t _p = 0; _p < (size_t)(NP); _p += 2) {                                    \
            const gfv _m = GFNI_TAIL_MAT4((BLK), (BLK), _p);                             \
            (t)[_p][GFNI_ZPP(pc) - 1] ^= gfv_affine(_zl2, _m);                               \
        }                                                                                    \
    } while (0)
/* the same from two plain last registers (each [e, 0, 0, 0 | o, 0, 0, 0]) */
#define GFNI_TAIL_QUAD(t, za, zc, BA, BB, NP)                                                \
    do {                                                                                     \
        const gfv _zl = _mm512_mask_blend_epi64(0xCC, _mm512_permutex_epi64((za), 0x00),     \
                                                _mm512_permutex_epi64((zc), 0x00));          \
        for (size_t _p = 0; _p < (size_t)(NP); _p += 2) {                                    \
            (t)[_p][GFNI_ZPP(pc) - 1] ^= gfv_affine(_zl, GFNI_TAIL_MAT4((BA), (BB), _p));    \
        }                                                                                    \
    } while (0)
/* t[p][last] and t[p+1][last] from the shared t[p][last] = [p, p+1, p, p+1 | ...]:
   the plain layout is [out, 0, 0, 0 | out', 0, 0, 0], as gfni_spread2 produces it */
#define GFNI_TAIL_UNPACK(t, NP)                                                              \
    for (size_t _p = 0; GFNI_LAST1 && _p < (size_t)(NP); _p += 2) {                          \
        const gfv _s = (t)[_p][GFNI_ZPP(pc) - 1];                                            \
        const gfv _g = _mm512_xor_si512(_s, _mm512_permutex_epi64(_s, 0x4E));                \
        (t)[_p][GFNI_ZPP(pc) - 1] = _mm512_maskz_mov_epi64(0x11, _g);                        \
        if (_p + 1 < (size_t)(NP))                                                           \
            (t)[_p + 1][GFNI_ZPP(pc) - 1] = _mm512_maskz_permutex_epi64(0x11, _g, 0x01);     \
    }
#else
/* ymm: only the 1-limb piece (MAYO_1, MAYO_5). Its spread pair (a, b) fills only qwords 0
   and 2, so (c, d) goes in qwords 1 and 3, with the matrix [M(a,b), M(c,d)] twice: four
   sources per affine instead of two. The pair steps leave qwords 1 and 3 zero, so the two
   mix freely; slots stored from such an accumulator hold the sum split over qwords 0 ^ 1
   and 2 ^ 3 (the consumers' affines are linear), and the despread folds it. */
#define GFNI_QUAD(pc) (GFNI_PC_LIMBS(pc) == 1)
#define GFNI_LAST1 1    /* a ymm quad step is always the packed 1-limb piece */
#define GFNI_TAILPACK 0
#define GFNI_TAIL_PAIR(t, zl, BLK, NP) ((void) 0)
#define GFNI_TAIL_QUAD(t, za, zc, BA, BB, NP) ((void) 0)
#define GFNI_TAIL_UNPACK(t, NP) ((void) 0)
#define GFNI_FOLD(z) _mm256_xor_si256((z), _mm256_srli_si256((z), 8))
static inline gfv gfni_bcast_limb(size_t pc, const uint64_t *p) {
    return _mm256_broadcastq_epi64(_mm_loadl_epi64((const __m128i *) (p + GFNI_PCL * pc)));
}
/* qwords [z0(a,b), z0(c,d), z1(a,b), z1(c,d)], z0/z1 as gfni_spread2 */
static inline gfv gfni_spread4(size_t pc, const uint64_t *a, const uint64_t *b, const uint64_t *c, const uint64_t *d) {
    const gfv x = _mm256_blend_epi32(gfni_bcast_limb(pc, a), gfni_bcast_limb(pc, c), 0xCC);
    const gfv y = _mm256_blend_epi32(gfni_bcast_limb(pc, b), gfni_bcast_limb(pc, d), 0xCC);
    return gfv_sel(gfv_set1_8(0x0f), _mm256_srlv_epi64(x, _mm256_setr_epi64x(0, 0, 4, 4)),
                   _mm256_sllv_epi64(y, _mm256_setr_epi64x(4, 4, 0, 0)));
}
#define GFNI_QUAD_STEP(t, a, b, c, d, BA, BB, NP)                                         \
    do {                                                                                     \
        const gfv _z = gfni_spread4(pc, (a), (b), (c), (d));                                 \
        for (size_t _p = 0; _p < (size_t)(NP); _p++) {                                       \
            gfv _m = _mm256_blend_epi32(gfv_set1_64((BA)[_p]), gfv_set1_64((BB)[_p]), 0xCC);                   \
            __asm__("" : "+v"(_m));                                                          \
            (t)[_p][0] ^= gfv_affine(_z, _m);                                                \
        }                                                                                    \
    } while (0)
#endif
static inline void gfni_spread2(size_t pc, const gfv x[GFNI_RP], const gfv y[GFNI_RP], gfv z[GFNI_ZPMAX]) {
    const gfv lo = gfv_set1_8(0x0f);
    for (size_t _n = 0; _n < GFNI_RP; _n++) {
        gfv _x = x[_n], _y = y[_n];
        gfv _x_shifted = gfv_nib_down(_x), _y_shifted = gfv_nib_up(_y);
        // select low nibble from x and high nibbles from y_shifted
        gfv z0 = gfv_sel(lo, _x, _y_shifted);
        // select low nibble from x_shifted and high nibbles from y
        gfv z1 = gfv_sel(lo, _x_shifted, _y);

        if (2 * _n + 1 < GFNI_ZPP(pc)) {
            z[2 * _n] = z0;
            z[2 * _n + 1] = z1;
        } else {
#if defined(MAYO_GFNI_ZMM) && GFNI_HALVES
            /* x and y fill both halves: z0 in the low half, z1 in the high one directly */
            (void) z0; (void) z1;
            z[2 * _n] = gfv_sel(lo, gfv_nib_down_hi(_x), gfv_nib_up_lo(_y));
#else
            z[2 * _n] = gfv_cat_lo(z0, z1);
#endif
        }
    }
}

/* one m-vec against a zero partner: in the low nibbles (HI = 0) or the high ones (HI = 1) */
static inline void gfni_spread1(size_t pc, const uint64_t *p, int hi, gfv z[GFNI_ZPMAX]) {
    gfv in[GFNI_RP], zero[GFNI_RP];
    gfni_load_mvec(pc, p, in);
    for (size_t q = 0; q < GFNI_RP; q++) zero[q] = gfv_zero();
    if (hi) gfni_spread2(pc, zero, in, z);
    else gfni_spread2(pc, in, zero, z);
}

static inline void gfni_despread2(size_t pc, const gfv z[GFNI_ZPMAX], gfv x[GFNI_RP], gfv y[GFNI_RP]) {
    const gfv lo = gfv_set1_8(0x0f);
    for (size_t _n = 0; _n < GFNI_RP; _n++) {
        gfv z0 = z[2 * _n], z1;
        if (2 * _n + 1 < GFNI_ZPP(pc)) {
            z1 = z[2 * _n + 1];
        } else {
            if (GFNI_QUAD(pc)) z0 = GFNI_FOLD(z0);
            z1 = gfv_dup_hi(z0);
        }
        gfv z0_shifted = gfv_nib_down(z0), z1_shifted = gfv_nib_up(z1);
        // select low nibble from x and high nibbles from y_shifted
        x[_n] = gfv_sel(lo, z0, z1_shifted);
        // select low nibble from x_shifted and high nibbles from y
        y[_n] = gfv_sel(lo, z0_shifted, z1);
    }
}
/* acc[0] ^= x and, if both, the m-vec after it ^= y */
static inline void gfni_xor_store_two(size_t pc, uint64_t *p, const gfv x[GFNI_RP], const gfv y[GFNI_RP], int both) {
#ifdef GFNI_PAIR_STORE_FITS
    if (both) {
        gfni_xor_store_pair(pc, p, x, y);
        return;
    }
#endif
    gfni_xor_store_mvec(pc, p, x);
    if (both) gfni_xor_store_mvec(pc, p + M_VEC_LIMBS_MAX, y);
}
// PS keeps the interleaved pairs, so nothing de-spreads until SPS: one pair slot is
// GFNI_SLOT limbs, KP of them to a row
#define MAYO_VKTMP_N (V_MAX * KP * GFNI_SLOT)
/* as KP, pairing along O_MAX */
#define OP ((O_MAX + 1) / 2)
#if defined(MAYO_GFNI_ZMM) && GFNI_TAILPACK && (KP % 2 || OP % 2)
#error "tail packing reads output pairs two at a time"
#endif
/* GFNI_BLK_TABLE for reduction-major tables, TABS[col * NOUT + out] */
#define GFNI_BLK_TABLE_RM(B, TABS, NOUT, NCOLS, NP)                                          \
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
                gfni_blk2(_a, _b, _c, _d);                                                   \
        }                                                                                    \
    }
#define PS2X2_PAIR(base, row, pair) ((base) + ((size_t)(row) * KP + (size_t)(pair)) * GFNI_SLOT)

/* lambda into column 0 of P1 * V^t. Column 0 is the even half of pair 0 here, so
   spread it against zero and xor the pair slot, the same full-width slot the
   producers store. v m-vector XORs, no multiplications. */
#define MAYO_HAVE_VKT_ADD_LAMBDA 1
/* for each row r < ROWS: the slot SLOT(r) ^= spread(lambda[r], 0) */
#define GFNI_ADD_EVEN_ROWS(ROWS, SLOT, LAMBDA)                                               \
    GFNI_FOR_PIECES(pc) {                                                                    \
        for (int r = 0; r < (ROWS); r++) {                                                   \
            uint64_t *slot = SLOT(r);                                                        \
            gfv z[GFNI_ZPMAX];                                                               \
            gfni_spread1(pc, (LAMBDA) + (size_t) r * M_VEC_LIMBS_MAX, 0, z);                 \
            for (size_t w = 0; w < GFNI_ZPP(pc); w++) {                                      \
                gfv_storeu(slot + GFNI_ZOFF(pc, w), gfv_xor(gfv_loadu(slot + GFNI_ZOFF(pc, w)), z[w])); \
            }                                                                                \
        }                                                                                    \
    }
static inline void m_vkt_add_lambda(uint64_t *VKtmp, const uint64_t *lambda, int v) {
#define VKT_SLOT(r) PS2X2_PAIR(VKtmp, r, 0)
    GFNI_ADD_EVEN_ROWS(v, VKT_SLOT, lambda);
#undef VKT_SLOT
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
static inline void mayo_gfni_V_multabs(const unsigned char *V, uint64_t *BV) {
    uint64_t V_tabs[K_MAX * V_MAX];
    mayo_gfni_tabs(V, V_tabs, (size_t) K_MAX * V_MAX);
    GFNI_BLK_TABLE(BV, V_tabs, V_TAB_STRIDE, V_MAX);
    mayo_secure_clear(V_tabs, sizeof(V_tabs));
}
// one 2x2 accumulate over both halves of an interleaved pair, NP output pairs. The empty
// asm keeps the matrix in a register: clang 14..19 miscompiles the folded m64bcst operand
// in unrolled loops, stepping its displacement by 64 bytes instead of 8.
#define GFNI_MACC2(t, z, blk, NP)                                                            \
    do {                                                                                     \
        GFNI_PRAGMA_(GCC unroll 16)                                                          \
        for (size_t _p = 0; _p < (size_t)(NP); _p++) {                                       \
            gfv _m = gfv_set1_64((blk)[_p]);                                                 \
            __asm__("" : "+v"(_m));                                                          \
            for (size_t _q = 0; _q + GFNI_TAILPACK < GFNI_ZPP(pc); _q++) {                    \
                (t)[_p][_q] ^= gfv_affine((z)[_q], _m);                                      \
            }                                                                                \
        }                                                                                    \
        if (GFNI_TAILPACK) GFNI_TAIL_PAIR(t, (z)[GFNI_ZPP(pc) - 1], (blk), NP);              \
    } while (0)
/* acc[l][j] (= or ^=) sum_r S[l][r] * PS[r][j] over ROWS contraction rows, for the
   products whose m-vec operand is stored as interleaved pairs (S^t * PS, V * Pv, U^t * O):
   those pairs run along the free index j, so a scalar fits only a 1x1 block there. Here
   column j's nibbles of two rows r, r+1 are put side by side instead, so one 2x2 block of
   S (outputs l, l+1 by rows r, r+1: B[rp * NP + lp], the table the P-side products use)
   does four products per affine. The accumulators are output pairs in the spread layout
   and despread as usual. Output l of column j goes to acc[l * SL + j * SJ]. */
#define GFNI_SPS_CH (GFNI_ACCS / GFNI_ZPMAX)   /* even: tail packing pairs output pairs */
/* column j's nibbles of rows 2rp, 2rp+1 side by side (column 2cp is the low nibbles of
   the stored pair, 2cp + 1 the high ones); a row past ROWS is zero */
#define GFNI_SPS_OPD(w, rp, ROWS)                                                                \
    do {                                                                                     \
        const uint64_t *_pa = PS2X2_PAIR(PS_, 2 * (rp), cp), *_pb = PS2X2_PAIR(PS_, 2 * (rp) + 1, cp); \
        for (size_t _q = 0; _q < GFNI_ZPP(pc); _q++) {                                       \
            const gfv _a = gfv_loadu(_pa + GFNI_ZOFF(pc, _q));                               \
            const gfv _b = 2 * (rp) + 1 < (size_t)(ROWS) ? gfv_loadu(_pb + GFNI_ZOFF(pc, _q)) : gfv_zero(); \
            (w)[_q] = (j & 1) ? gfv_sel(lo, gfv_nib_down(_a), _b) : gfv_sel(lo, _a, gfv_nib_up(_b)); \
        }                                                                                    \
    } while (0)
#define GFNI_SPS(PS, ROWS, B, NOUT, ACC, SL, SJ, ACCUMULATE)                                  \
    do {                                                                                     \
    const uint64_t *const PS_ = (PS), *const B_ = (B); uint64_t *const acc_ = (ACC);          \
    const gfv lo = gfv_set1_8(0x0f);                                                         \
    const size_t np = ((NOUT) + 1) / 2;                                                      \
    GFNI_FOR_PIECES(pc) {                                                                    \
        for (size_t j = 0; j < K_MAX; j++) {                                                 \
            const size_t cp = j / 2;                                                         \
            for (size_t lp0 = 0; lp0 < np; lp0 += GFNI_SPS_CH) {                             \
                const size_t lpn = np - lp0 < (size_t) GFNI_SPS_CH ? np - lp0 : (size_t) GFNI_SPS_CH; \
                gfv t[GFNI_SPS_CH][GFNI_ZPMAX];                                              \
                for (size_t l = 0; l < lpn; l++) {                                           \
                    for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[l][q] = gfv_zero();          \
                }                                                                            \
                size_t rp = 0;                                                               \
                for (; 2 * rp + 3 < (size_t)(ROWS); rp += 2) {   /* two row pairs per pass */  \
                    gfv w[GFNI_ZPMAX], w2[GFNI_ZPMAX];                                       \
                    GFNI_SPS_OPD(w, rp, ROWS); GFNI_SPS_OPD(w2, rp + 1, ROWS);                           \
                    for (size_t l = 0; l < lpn; l++) {                                       \
                        gfv m = gfv_set1_64(B_[rp * np + lp0 + l]), m2 = gfv_set1_64(B_[(rp + 1) * np + lp0 + l]); \
                        __asm__("" : "+v"(m), "+v"(m2));                                     \
                        for (size_t q = 0; q + GFNI_TAILPACK < GFNI_ZPP(pc); q++)            \
                            t[l][q] = gfv_xor3(t[l][q], gfv_affine(w[q], m), gfv_affine(w2[q], m2)); \
                    }                                                                        \
                    if (GFNI_TAILPACK) GFNI_TAIL_QUAD(t, w[GFNI_ZPP(pc) - 1], w2[GFNI_ZPP(pc) - 1], \
                                                      B_ + rp * np + lp0, B_ + (rp + 1) * np + lp0, lpn); \
                }                                                                            \
                for (; 2 * rp < (size_t)(ROWS); rp++) {                                      \
                    gfv w[GFNI_ZPMAX];                                                       \
                    GFNI_SPS_OPD(w, rp, ROWS);                                                 \
                    GFNI_MACC2(t, w, B_ + rp * np + lp0, lpn);                               \
                }                                                                            \
                GFNI_TAIL_UNPACK(t, lpn);                                                    \
                for (size_t l = 0; l < lpn; l++) {                                           \
                    const size_t o = 2 * (lp0 + l);                                          \
                    uint64_t *d0 = acc_ + (o * (SL) + j * (SJ)) * M_VEC_LIMBS_MAX, *d1 = d0 + (SL) * M_VEC_LIMBS_MAX; \
                    gfv x[GFNI_RP], y[GFNI_RP];                                              \
                    gfni_despread2(pc, t[l], x, y);                                          \
                    if (ACCUMULATE) {                                                        \
                        gfni_xor_store_mvec(pc, d0, x);                                      \
                        if (o + 1 < (size_t)(NOUT)) gfni_xor_store_mvec(pc, d1, y);          \
                    } else {                                                                 \
                        gfni_store_mvec(pc, d0, x);                                          \
                        if (o + 1 < (size_t)(NOUT)) gfni_store_mvec(pc, d1, y);              \
                    }                                                                        \
                }                                                                            \
            }                                                                                \
        }                                                                                    \
    }                                                                                        \
    } while (0)




// P1 * O -> P1: v x v (upper triangular), O: v x o
static
inline void P1_times_O(const uint64_t *P1, const uint64_t *O_tabs, uint64_t *acc) {
    uint64_t BO[((V_MAX + 1) / 2) * OP];
    GFNI_BLK_TABLE_RM(BO, O_tabs, O_MAX, V_MAX, OP);
    GFNI_FOR_PIECES(pc) {
        size_t cols_used = 0;
        for (size_t r = 0; r < V_MAX; r++) {
            gfv t[OP][GFNI_ZPMAX];
            for (size_t lp = 0; lp < OP; lp++) {
                for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[lp][q] = gfv_zero();
            }
            size_t c = r;
            if (c & 1) {            /* an odd start puts the first column in the pair's second half */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P1, cols_used), 1, z);
                cols_used++;
                GFNI_MACC2(t, z, BO + (c / 2) * OP, OP);
                c++;
            }
            for (; GFNI_QUAD(pc) && c + 3 < V_MAX; c += 4, cols_used += 4) {
                GFNI_QUAD_STEP(t, PINV(P1, cols_used), PINV(P1, cols_used + 1), PINV(P1, cols_used + 2),
                               PINV(P1, cols_used + 3), BO + (c / 2) * OP, BO + (c / 2 + 1) * OP, OP);
            }
            for (; c + 1 < V_MAX; c += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, PINV(P1, cols_used), a);
                gfni_load_mvec(pc, PINV(P1, cols_used + 1), b);
                cols_used += 2;
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, BO + (c / 2) * OP, OP);
            }
            if (c < V_MAX) {        /* odd V_MAX */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P1, cols_used), 0, z);
                cols_used++;
                GFNI_MACC2(t, z, BO + (c / 2) * OP, OP);
            }
            GFNI_TAIL_UNPACK(t, OP);
            for (size_t lp = 0; lp < OP; lp++) {
                gfv x[GFNI_RP], y[GFNI_RP];
                gfni_despread2(pc, t[lp], x, y);
                gfni_xor_store_two(pc, acc + (r * O_MAX + 2 * lp) * M_VEC_LIMBS_MAX, x, y, 2 * lp + 1 < O_MAX);
            }
        }
    }
    mayo_secure_clear(BO, sizeof(BO));
}

static
inline void Ot_times_P1O_P2(const uint64_t *P1O_P2, const uint64_t *O_tabs, uint64_t *acc) {
    uint64_t BO[((V_MAX + 1) / 2) * OP];
    GFNI_BLK_TABLE_RM(BO, O_tabs, O_MAX, V_MAX, OP);
    GFNI_FOR_PIECES(pc) {
        for (size_t c = 0; c < O_MAX; c++) {
            gfv t[OP][GFNI_ZPMAX];
            for (size_t lp = 0; lp < OP; lp++) {
                for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[lp][q] = gfv_zero();
            }
            size_t r = 0;
#define P1OP2(r) (P1O_P2 + ((r) * O_MAX + c) * M_VEC_LIMBS_MAX)
            for (; GFNI_QUAD(pc) && r + 3 < V_MAX; r += 4) {
                GFNI_QUAD_STEP(t, P1OP2(r), P1OP2(r + 1), P1OP2(r + 2), P1OP2(r + 3),
                               BO + (r / 2) * OP, BO + (r / 2 + 1) * OP, OP);
            }
#undef P1OP2
            for (; r + 1 < V_MAX; r += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, P1O_P2 + (r * O_MAX + c) * M_VEC_LIMBS_MAX, a);
                gfni_load_mvec(pc, P1O_P2 + ((r + 1) * O_MAX + c) * M_VEC_LIMBS_MAX, b);
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, BO + (r / 2) * OP, OP);
            }
            if (r < V_MAX) {        /* odd V_MAX */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, P1O_P2 + (r * O_MAX + c) * M_VEC_LIMBS_MAX, 0, z);
                GFNI_MACC2(t, z, BO + (r / 2) * OP, OP);
            }
            GFNI_TAIL_UNPACK(t, OP);
            for (size_t lp = 0; lp < OP; lp++) {
                gfv x[GFNI_RP], y[GFNI_RP];
                gfni_despread2(pc, t[lp], x, y);
                gfni_xor_store_mvec(pc, acc + ((2 * lp) * O_MAX + c) * M_VEC_LIMBS_MAX, x);
                if (2 * lp + 1 < O_MAX) {
                    gfni_xor_store_mvec(pc, acc + ((2 * lp + 1) * O_MAX + c) * M_VEC_LIMBS_MAX, y);
                }
            }
        }
    }
    mayo_secure_clear(BO, sizeof(BO));
}

static
inline void P1P1t_times_O(const mayo_params_t *p, const uint64_t *P1, const unsigned char *O, uint64_t *acc) {
    (void) p;
    uint64_t O_tabs[MAYO_O_MULTABS_N];
    mayo_O_multabs(O, O_tabs);
    uint64_t BO[((V_MAX + 1) / 2) * OP];
    GFNI_BLK_TABLE_RM(BO, O_tabs, O_MAX, V_MAX, OP);
    GFNI_FOR_PIECES(pc) {
        size_t diag = 0;            /* packed index of (r, r) */
        for (size_t r = 0; r < V_MAX; r++) {
            gfv t[OP][GFNI_ZPMAX];
            for (size_t lp = 0; lp < OP; lp++) {
                for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[lp][q] = gfv_zero();
            }
            /* columns left of r come down column r, those right of r along row r; column r
               drops out, P1 + P1t having a zero diagonal. The split keeps its pair out of
               both runs. */
            const size_t rp = r & ~(size_t)1;
            size_t pos = r;         /* packed index of (c, r) */
            size_t c = 0;
            for (; GFNI_QUAD(pc) && c + 4 <= rp; c += 4) {
                const size_t p1 = pos + (V_MAX - c - 1), p2 = p1 + (V_MAX - c - 2), p3 = p2 + (V_MAX - c - 3);
                GFNI_QUAD_STEP(t, PINV(P1, pos), PINV(P1, p1), PINV(P1, p2), PINV(P1, p3),
                               BO + (c / 2) * OP, BO + (c / 2 + 1) * OP, OP);
                pos = p3 + (V_MAX - c - 4);
            }
            for (; c < rp; c += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, PINV(P1, pos), a);
                pos += (V_MAX - c - 1);
                gfni_load_mvec(pc, PINV(P1, pos), b);
                pos += (V_MAX - c - 2);
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, BO + (c / 2) * OP, OP);
            }
            {                       /* the pair holding the diagonal: one column, or none */
                gfv z[GFNI_ZPMAX];
                if (r & 1) {
                    gfni_spread1(pc, PINV(P1, pos), 0, z);   /* (r - 1, r) */
                    GFNI_MACC2(t, z, BO + (rp / 2) * OP, OP);
                } else if (r + 1 < V_MAX) {
                    gfni_spread1(pc, PINV(P1, diag + 1), 1, z);
                    GFNI_MACC2(t, z, BO + (rp / 2) * OP, OP);
                }
                c = rp + 2;
            }
            for (; GFNI_QUAD(pc) && c + 3 < V_MAX; c += 4) {
                GFNI_QUAD_STEP(t, PINV(P1, diag + (c - r)), PINV(P1, diag + (c + 1 - r)), PINV(P1, diag + (c + 2 - r)),
                               PINV(P1, diag + (c + 3 - r)), BO + (c / 2) * OP, BO + (c / 2 + 1) * OP, OP);
            }
            for (; c + 1 < V_MAX; c += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, PINV(P1, diag + (c - r)), a);
                gfni_load_mvec(pc, PINV(P1, diag + (c + 1 - r)), b);
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, BO + (c / 2) * OP, OP);
            }
            if (c < V_MAX) {        /* odd V_MAX */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P1, diag + (c - r)), 0, z);
                GFNI_MACC2(t, z, BO + (c / 2) * OP, OP);
            }
            GFNI_TAIL_UNPACK(t, OP);
            for (size_t lp = 0; lp < OP; lp++) {
                gfv x[GFNI_RP], y[GFNI_RP];
                gfni_despread2(pc, t[lp], x, y);
                gfni_xor_store_two(pc, acc + (r * O_MAX + 2 * lp) * M_VEC_LIMBS_MAX, x, y, 2 * lp + 1 < O_MAX);
            }
            diag += (V_MAX - r);
        }
    }
    mayo_secure_clear(O_tabs, sizeof(O_tabs));
    mayo_secure_clear(BO, sizeof(BO));
}

static
inline void Vt_times_L(const uint64_t *L, const uint64_t *BV, uint64_t *acc) {
    /* L is in the plain layout, so the pair is spread here, not loaded interleaved */
    GFNI_FOR_PIECES(pc) {
        for (size_t c = 0; c < O_MAX; c++) {
            gfv t[KP][GFNI_ZPMAX];
            for (size_t lp = 0; lp < KP; lp++) {
                for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[lp][q] = gfv_zero();
            }
            size_t r = 0;
#define LRC(r) (L + ((r) * O_MAX + c) * M_VEC_LIMBS_MAX)
            for (; GFNI_QUAD(pc) && r + 3 < V_MAX; r += 4) {
                GFNI_QUAD_STEP(t, LRC(r), LRC(r + 1), LRC(r + 2), LRC(r + 3),
                               BV + (r / 2) * KP, BV + (r / 2 + 1) * KP, KP);
            }
#undef LRC
            for (; r + 1 < V_MAX; r += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, L + (r * O_MAX + c) * M_VEC_LIMBS_MAX, a);
                gfni_load_mvec(pc, L + ((r + 1) * O_MAX + c) * M_VEC_LIMBS_MAX, b);
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, BV + (r / 2) * KP, KP);
            }
            if (r < V_MAX) {            /* odd V_MAX: the last row pairs with zero */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, L + (r * O_MAX + c) * M_VEC_LIMBS_MAX, 0, z);
                GFNI_MACC2(t, z, BV + (r / 2) * KP, KP);
            }
            GFNI_TAIL_UNPACK(t, KP);
            for (size_t lp = 0; lp < KP; lp++) {
                gfv x[GFNI_RP], y[GFNI_RP];
                gfni_despread2(pc, t[lp], x, y);
                gfni_xor_store_mvec(pc, acc + ((2 * lp) * O_MAX + c) * M_VEC_LIMBS_MAX, x);
                if (2 * lp + 1 < K_MAX) {
                    gfni_xor_store_mvec(pc, acc + ((2 * lp + 1) * O_MAX + c) * M_VEC_LIMBS_MAX, y);
                }
            }
        }
    }
}

static
inline void Vt_times_Pv(const uint64_t *Pv, const uint64_t *BV, uint64_t *acc) {
    GFNI_SPS(Pv, V_MAX, BV, K_MAX, acc, K_MAX, 1, 0);
}

static
inline void P1_times_Vt(const uint64_t *P1, const uint64_t *BV, uint64_t *acc) {
    GFNI_FOR_PIECES(pc) {
        size_t cols_used = 0;
        for (size_t r = 0; r < V_MAX; r++) {
            gfv t[KP][GFNI_ZPMAX];
            for (size_t lp = 0; lp < KP; lp++) {
                for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[lp][q] = gfv_zero();
            }
            size_t c = r;
            if (c & 1) {
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P1, cols_used), 1, z);
                cols_used++;
                GFNI_MACC2(t, z, BV + (c / 2) * KP, KP);
                c++;
            }
            for (; GFNI_QUAD(pc) && c + 3 < V_MAX; c += 4, cols_used += 4) {
                GFNI_QUAD_STEP(t, PINV(P1, cols_used), PINV(P1, cols_used + 1), PINV(P1, cols_used + 2),
                               PINV(P1, cols_used + 3), BV + (c / 2) * KP, BV + (c / 2 + 1) * KP, KP);
            }
            for (; c + 1 < V_MAX; c += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, PINV(P1, cols_used), a);
                gfni_load_mvec(pc, PINV(P1, cols_used + 1), b);
                cols_used += 2;
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, BV + (c / 2) * KP, KP);
            }
            if (c < V_MAX) {        /* odd V_MAX: a last column with no partner */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P1, cols_used), 0, z);
                cols_used++;
                GFNI_MACC2(t, z, BV + (c / 2) * KP, KP);
            }
            GFNI_TAIL_UNPACK(t, KP);
            for (size_t lp = 0; lp < KP; lp++) {
                uint64_t *q = PS2X2_PAIR(acc, r, lp);
                for (size_t w = 0; w < GFNI_ZPP(pc); w++) gfv_storeu((q + GFNI_ZOFF(pc, w)), t[lp][w]);
            }
        }
    }
}

// acc += P1^t * V^t, where P1 is stored as an upper-triangular matrix.
static
inline void P1t_times_Vt(const uint64_t *P1, const uint64_t *BV, uint64_t *acc) {
    GFNI_FOR_PIECES(pc) {
        for (size_t r = 0; r < V_MAX; r++) {
            gfv t[KP][GFNI_ZPMAX];
            for (size_t lp = 0; lp < KP; lp++) {
                for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[lp][q] = gfv_zero();
            }
            /* the walk lands exactly on the diagonal element, so one loop covers c = 0..r */
            size_t pos = r, c = 0;
            for (; GFNI_QUAD(pc) && c + 3 <= r; c += 4) {
                const size_t p1 = pos + (V_MAX - c - 1), p2 = p1 + (V_MAX - c - 2), p3 = p2 + (V_MAX - c - 3);
                GFNI_QUAD_STEP(t, PINV(P1, pos), PINV(P1, p1), PINV(P1, p2), PINV(P1, p3),
                               BV + (c / 2) * KP, BV + (c / 2 + 1) * KP, KP);
                pos = p3 + (V_MAX - c - 4);
            }
            for (; c + 1 <= r; c += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, PINV(P1, pos), a);
                pos += (V_MAX - c - 1);
                gfni_load_mvec(pc, PINV(P1, pos), b);
                pos += (V_MAX - c - 2);
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, BV + (c / 2) * KP, KP);
            }
            if (c <= r) {
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P1, pos), 0, z);
                GFNI_MACC2(t, z, BV + (c / 2) * KP, KP);
            }
            GFNI_TAIL_UNPACK(t, KP);
            for (size_t lp = 0; lp < KP; lp++) {
                uint64_t *q = PS2X2_PAIR(acc, r, lp);
                for (size_t w = 0; w < GFNI_ZPP(pc); w++) {
                    gfv_storeu((q + GFNI_ZOFF(pc, w)),
                        gfv_loadu((q + GFNI_ZOFF(pc, w))) ^ t[lp][w]);
                }
            }
        }
    }
}

/* acc[c] += sum_r O[r][c] * lambda[r]: the one column of U^t * O that the
   materialized-L signing path has to supply itself. Same shape as Ut_times_O
   below with a single column pair and only its even half live, so lambda is
   read unspread and nothing has to be de-spread on the way out. */
#define MAYO_HAVE_LAMBDA_TIMES_O 1
static
inline void m_lambda_times_O(const uint64_t *lambda, const mayo_multab_t *O_tabs,
                                    uint64_t *acc) {
    GFNI_FOR_PIECES(pc) {
        for (size_t l0 = 0; l0 < O_MAX; l0 += OCH) {
            const size_t ln = (O_MAX - l0 < (size_t) OCH) ? O_MAX - l0 : (size_t) OCH;
            gfv t[OCH][GFNI_RP];
            for (size_t li = 0; li < ln; li++) {
                for (size_t q = 0; q < GFNI_RP; q++) {
                    t[li][q] = gfv_zero();
                }
            }
            for (size_t r = 0; r < V_MAX; r++) {
                gfv z[GFNI_RP];
                gfni_load_mvec(pc, lambda + r * M_VEC_LIMBS_MAX, z);
                for (size_t li = 0; li < ln; li++) {
                    gfv m = gfv_set1_64(O_tabs[O_MAX * r + l0 + li]);
                    __asm__("" : "+v"(m));
                    for (size_t q = 0; q < GFNI_RP; q++) {
                        t[li][q] ^= gfv_affine(z[q], m);
                    }
                }
            }
            for (size_t li = 0; li < ln; li++) {
                gfni_xor_store_mvec(pc, acc + (l0 + li) * M_VEC_LIMBS_MAX, t[li]);
            }
        }
    }
}

// acc += U^t * O, where U is a v-by-k matrix of m-vectors.
static
inline void Ut_times_O(const uint64_t *U, const uint64_t *O_tabs, uint64_t *acc) {
    uint64_t BO[((V_MAX + 1) / 2) * OP];
    GFNI_BLK_TABLE_RM(BO, O_tabs, O_MAX, V_MAX, OP);
    GFNI_SPS(U, V_MAX, BO, O_MAX, acc, 1, O_MAX, 1);
    mayo_secure_clear(BO, sizeof(BO));
}


// PS rows 0..v-1: P1 * S1 + P2 * S2
static
inline void P1_times_S1_plus_P2_times_S2(const uint64_t *P1, const uint64_t *P2, const uint64_t *B1,
                                         const uint64_t *B2, uint64_t *acc) {
    /* One affine per input pair and output pair; an odd start spreads (0, x) so the
       pair entry still applies. */
    GFNI_FOR_PIECES(pc) {

        size_t P1_vec = 0;
        for (size_t r = 0; r < V_MAX; r++) {
            gfv t[KP][GFNI_ZPMAX];
            for (size_t lp = 0; lp < KP; lp++) {
                for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[lp][q] = gfv_zero();
            }

            // P1 * S1
            size_t c = r;
            if (c & 1) {
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P1, P1_vec), 1, z);
                P1_vec++;
                GFNI_MACC2(t, z, B1 + (c / 2) * KP, KP);
                c++;
            }
            for (; GFNI_QUAD(pc) && c + 3 < V_MAX; c += 4, P1_vec += 4) {
                GFNI_QUAD_STEP(t, PINV(P1, P1_vec), PINV(P1, P1_vec + 1), PINV(P1, P1_vec + 2),
                               PINV(P1, P1_vec + 3), B1 + (c / 2) * KP, B1 + (c / 2 + 1) * KP, KP);
            }
            for (; c + 1 < V_MAX; c += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, PINV(P1, P1_vec), a);
                gfni_load_mvec(pc, PINV(P1, P1_vec + 1), b);
                P1_vec += 2;
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, B1 + (c / 2) * KP, KP);
            }
            if (c < V_MAX) {        /* odd V_MAX: a last column with no partner */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P1, P1_vec), 0, z);
                P1_vec++;
                GFNI_MACC2(t, z, B1 + (c / 2) * KP, KP);
            }

            // P2 * S2
            size_t cc = 0;
            for (; cc + 1 < O_MAX; cc += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, PINV(P2, r * O_MAX + cc), a);
                gfni_load_mvec(pc, PINV(P2, r * O_MAX + cc + 1), b);
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, B2 + (cc / 2) * KP, KP);
            }
            if (cc < O_MAX) {       /* odd O_MAX */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, PINV(P2, r * O_MAX + cc), 0, z);
                GFNI_MACC2(t, z, B2 + (cc / 2) * KP, KP);
            }

            GFNI_TAIL_UNPACK(t, KP);
            for (size_t lp = 0; lp < KP; lp++) {
                uint64_t *q = PS2X2_PAIR(acc, r, lp);
                for (size_t w = 0; w < GFNI_ZPP(pc); w++) gfv_storeu((q + GFNI_ZOFF(pc, w)), t[lp][w]);
            }
        }
    }
}

// P3*S2 -> P3: o x o, S2: o x k // P3 upper triangular
static
inline void P3_times_S2(const uint64_t *P3, const uint64_t *B, uint64_t *acc) {
    GFNI_FOR_PIECES(pc) {
        size_t cols_used = 0;
        for (size_t r = 0; r < O_MAX; r++) {
            gfv t[KP][GFNI_ZPMAX];
            for (size_t lp = 0; lp < KP; lp++) {
                for (size_t q = 0; q < GFNI_ZPP(pc); q++) t[lp][q] = gfv_zero();
            }
            size_t c = r;
            if (c & 1) {
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, P3 + cols_used, 1, z);
                cols_used += M_VEC_LIMBS_MAX;
                GFNI_MACC2(t, z, B + (c / 2) * KP, KP);
                c++;
            }
            for (; GFNI_QUAD(pc) && c + 3 < O_MAX; c += 4, cols_used += 4 * M_VEC_LIMBS_MAX) {
                GFNI_QUAD_STEP(t, P3 + cols_used, P3 + cols_used + M_VEC_LIMBS_MAX, P3 + cols_used + 2 * M_VEC_LIMBS_MAX,
                               P3 + cols_used + 3 * M_VEC_LIMBS_MAX, B + (c / 2) * KP, B + (c / 2 + 1) * KP, KP);
            }
            for (; c + 1 < O_MAX; c += 2) {
                gfv a[GFNI_RP], b[GFNI_RP], z[GFNI_ZPMAX];
                gfni_load_mvec(pc, P3 + cols_used, a);
                cols_used += M_VEC_LIMBS_MAX;
                gfni_load_mvec(pc, P3 + cols_used, b);
                cols_used += M_VEC_LIMBS_MAX;
                gfni_spread2(pc, a, b, z);
                GFNI_MACC2(t, z, B + (c / 2) * KP, KP);
            }
            if (c < O_MAX) {        /* odd O_MAX */
                gfv z[GFNI_ZPMAX];
                gfni_spread1(pc, P3 + cols_used, 0, z);
                cols_used += M_VEC_LIMBS_MAX;
                GFNI_MACC2(t, z, B + (c / 2) * KP, KP);
            }
            GFNI_TAIL_UNPACK(t, KP);
            for (size_t lp = 0; lp < KP; lp++) {
                uint64_t *q = PS2X2_PAIR(acc, r, lp);
                for (size_t w = 0; w < GFNI_ZPP(pc); w++) gfv_storeu((q + GFNI_ZOFF(pc, w)), t[lp][w]);
            }
        }
    }
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
    for (int r = 0; r < k; r++) {   /* split each row of S at the vinegar/oil boundary */
        memcpy(S1 + (size_t) r * v, S + (size_t) r * n, (size_t) v);
        memcpy(S2 + (size_t) r * o, S + (size_t) r * n + v, (size_t) o);
    }

    alignas(64) uint64_t PS[N_MAX * KP * GFNI_SLOT]; // interleaved pairs, fully written below
#define MAYO_PS_ROWS_OFF (V_MAX * KP * GFNI_SLOT)

    uint64_t S1_tabs[MAYO_S1_MULTABS_N];
    uint64_t S2_tabs[MAYO_S2_MULTABS_N];
    mayo_gfni_tabs(S1, S1_tabs, (size_t) K_MAX * V_MAX);
    mayo_gfni_tabs(S2, S2_tabs, (size_t) K_MAX * O_MAX);
    /* 2x2 blocks of S1 and S2, shared by P * S^t and S^t * PS */
    uint64_t B1[((V_MAX + 1) / 2) * KP], B2[((O_MAX + 1) / 2) * KP];
    GFNI_BLK_TABLE(B1, S1_tabs, V_TAB_STRIDE, V_MAX);
    GFNI_BLK_TABLE(B2, S2_tabs, S2_TAB_STRIDE, O_MAX);

    P1_times_S1_plus_P2_times_S2(P1, P2, B1, B2, PS);
    P3_times_S2(P3, B2, PS + MAYO_PS_ROWS_OFF); // upper triangular

    /* lambda occupies column 0 of PS, which in this layout is the even half of column
       pair 0; spread it into that half and xor it into the pair slot, the same full-width
       slot the producers above store. The S^t * PS pass then applies S[i][r] to it, so
       SPS[i][0] picks up Lambda(s_i) for n m-vector XORs and no multiplications. */
#define PS_SLOT(r) ((r) < v ? PS2X2_PAIR(PS, r, 0) : PS2X2_PAIR(PS + MAYO_PS_ROWS_OFF, (r) - v, 0))
    GFNI_ADD_EVEN_ROWS(n, PS_SLOT, lambda);
#undef PS_SLOT

    // S^T * PS = S1^t*PS1 + S2^t*PS2
    GFNI_SPS(PS, V_MAX, B1, K_MAX, SPS, K_MAX, 1, 0);
    GFNI_SPS(PS + MAYO_PS_ROWS_OFF, O_MAX, B2, K_MAX, SPS, K_MAX, 1, 1);
}


#ifdef MAYO_GFNI_ZMM
/* transpose_16x16_nibbles with the 16 rows in two zmm. Each butterfly stage pairs rows d
   apart (d = 1, 2, 4 within a register, 8 across): t = ((row >> w) ^ partner) & mask is
   formed on every lane against the permuted partner, then the row with bit d clear takes
   t << w and its partner takes t, moved back by the same permute. */
#define MAYO_HAVE_TRANSPOSE16 1
#define GFNI_T16_STAGE(z, w, mask, PERM, LOW)                                                \
    do {                                                                                     \
        const __m512i _t = _mm512_ternarylogic_epi64(_mm512_srli_epi64((z), (w)), PERM(z),   \
                                                     _mm512_set1_epi64((long long) (mask)), 0x28); \
        const __m512i _u = PERM(_t);    /* t of the low row, on the partner's lane */       \
        (z) = _mm512_xor_si512((z), _mm512_mask_blend_epi64((LOW), _u, _mm512_slli_epi64(_t, (w)))); \
    } while (0)
#define GFNI_T16_P1(z) _mm512_shuffle_epi32((z), (_MM_PERM_ENUM) 0x4E)
#define GFNI_T16_P2(z) _mm512_permutex_epi64((z), 0x4E)
#define GFNI_T16_P4(z) _mm512_shuffle_i64x2((z), (z), 0x4E)
static inline void transpose_16x16_nibbles(uint64_t *M) {
    __m512i z0 = _mm512_loadu_si512((const void *) M), z1 = _mm512_loadu_si512((const void *) (M + 8));
    GFNI_T16_STAGE(z0, 4, 0x0f0f0f0f0f0f0f0fULL, GFNI_T16_P1, 0x55);
    GFNI_T16_STAGE(z1, 4, 0x0f0f0f0f0f0f0f0fULL, GFNI_T16_P1, 0x55);
    GFNI_T16_STAGE(z0, 8, 0x00ff00ff00ff00ffULL, GFNI_T16_P2, 0x33);
    GFNI_T16_STAGE(z1, 8, 0x00ff00ff00ff00ffULL, GFNI_T16_P2, 0x33);
    GFNI_T16_STAGE(z0, 16, 0x0000ffff0000ffffULL, GFNI_T16_P4, 0x0F);
    GFNI_T16_STAGE(z1, 16, 0x0000ffff0000ffffULL, GFNI_T16_P4, 0x0F);
    const __m512i t = _mm512_ternarylogic_epi64(_mm512_srli_epi64(z0, 32), z1,
                                                _mm512_set1_epi64(0x00000000ffffffffLL), 0x28);
    _mm512_storeu_si512((void *) M, _mm512_xor_si512(z0, _mm512_slli_epi64(t, 32)));
    _mm512_storeu_si512((void *) (M + 8), _mm512_xor_si512(z1, t));
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

#endif

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

/* dst[i] ^= each nibble of src[i] times coeff, one affine multiply per register */
#define MAYO_HAVE_GF16_MULC_XOR 1
static inline void m_gf16_mulc_xor(uint64_t *dst, const uint64_t *src, int n, unsigned char coeff, const unsigned char *t16) {
    (void) t16;   /* the affine matrix folds away for a constant coeff */
    const gfv mat = gfv_set1_64(mayo_gfni_tab(coeff));
    const int batched = n & ~(GFNI_W - 1);
    int i = 0;
    for (; i < batched; i += GFNI_W) {
        gfv v = gfv_loadu(src + i);
        gfv_storeu(dst + i, gfv_xor(gfv_loadu(dst + i), gfv_affine(v, mat)));
    }
    if (i < n) {
#ifdef MAYO_GFNI_ZMM
        const __mmask8 km = (__mmask8) ((1u << (n - i)) - 1u);
        gfv v = _mm512_maskz_loadu_epi64(km, src + i);
        gfv d = _mm512_maskz_loadu_epi64(km, dst + i);
        _mm512_mask_storeu_epi64(dst + i, km, gfv_xor(d, gfv_affine(v, mat)));
#else
        const __m256i km = _mm256_cmpgt_epi64(_mm256_set1_epi64x(n - i), _mm256_setr_epi64x(0, 1, 2, 3));
        gfv v = _mm256_maskload_epi64((const long long *) (src + i), km);
        gfv d = _mm256_maskload_epi64((const long long *) (dst + i), km);
        _mm256_maskstore_epi64((long long *) (dst + i), km, d ^ gfv_affine(v, mat));
#endif
    }
}

/* see the fallback in arithmetic.h; used per use on ymm (avx512_gfni takes all uses of a
   block at once, m_compute_A_block below). Four columns are loaded whole and transposed
   in registers rather than gathered, then two, then one; the loads run past the last
   column into Mtmp's MAYO_MVEC_SLACK tail. Every index and the shift amount are public. */
#define MAYO_HAVE_COMPUTE_A_ACC 1
static inline void m_compute_A_acc(uint64_t *A, size_t A_width, int row0, int col,
                                   const uint64_t *M, int o, int mvl, int bits) {
    const __m128i shl = _mm_cvtsi32_si128(bits);
    const __m128i shr = _mm_cvtsi32_si128(64 - bits);
    int c = 0;
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

#ifdef MAYO_GFNI_ZMM
/* All uses of one o-by-mvl block M at once: transpose it once (Mt row l =
   limb l of the o m-vecs), then per use one read-modify-write per destination row and
   column chunk, limb l shifted up and limb l-1's carry in one ternlog. Whole chunks of
   8, 4, 2 columns with plain loads and stores (masked ones split badly on A's
   unaligned rows), then at most one column. */
#define MAYO_HAVE_COMPUTE_A_BLOCK 1
#define GFNI_MT_S (((O_MAX + 7) / 8) * 8)
/* one use on a W-column chunk at cc: V, LD/ST/XOR/SLL/SRL/TERN of that width */
#define GFNI_A_USE(V, LD, ST, XOR, SLL, SRL, TERN)                                           \
    do {                                                                                     \
        uint64_t *_d = A + (size_t) row0[u] * A_width + col[u] + cc;                         \
        V _prev = LD(Mt + cc);                                                               \
        ST(_d, XOR(LD(_d), SLL(_prev, shl)));                                                \
        if (bits[u] > 0) {                                                                   \
            for (int _l = 1; _l < mvl; _l++) {                                               \
                const V _cur = LD(Mt + (size_t) _l * GFNI_MT_S + cc);                        \
                uint64_t *_dl = _d + (size_t) _l * A_width;                                  \
                ST(_dl, TERN(LD(_dl), SLL(_cur, shl), SRL(_prev, shr)));                     \
                _prev = _cur;                                                                \
            }                                                                                \
            uint64_t *_dt = _d + (size_t) mvl * A_width;                                     \
            ST(_dt, XOR(LD(_dt), SRL(_prev, shr)));                                          \
        } else {                                                                             \
            for (int _l = 1; _l < mvl; _l++) {                                               \
                uint64_t *_dl = _d + (size_t) _l * A_width;                                  \
                ST(_dl, XOR(LD(_dl), LD(Mt + (size_t) _l * GFNI_MT_S + cc)));                \
            }                                                                                \
        }                                                                                    \
    } while (0)
#define GFNI_LD8(p) _mm512_loadu_si512((const void *) (p))
#define GFNI_ST8(p, x) _mm512_storeu_si512((void *) (p), (x))
#define GFNI_TERN8(a, b, c) _mm512_ternarylogic_epi64((a), (b), (c), 0x96)
#define GFNI_LD4(p) _mm256_loadu_si256((const __m256i *) (p))
#define GFNI_ST4(p, x) _mm256_storeu_si256((__m256i *) (p), (x))
#define GFNI_TERN4(a, b, c) _mm256_xor_si256((a), _mm256_xor_si256((b), (c)))
#define GFNI_LD2(p) _mm_loadu_si128((const __m128i *) (p))
#define GFNI_ST2(p, x) _mm_storeu_si128((__m128i *) (p), (x))
#define GFNI_TERN2(a, b, c) _mm_xor_si128((a), _mm_xor_si128((b), (c)))
static inline void m_compute_A_block(uint64_t *A, size_t A_width, const uint64_t *M, int o, int mvl,
                                     int nuse, const int *row0, const int *col, const int *bits) {
    alignas(64) uint64_t Mt[M_VEC_LIMBS_MAX * GFNI_MT_S];
    int c = 0;
    for (; c + 8 <= o; c += 8) {
        for (int l0 = 0; l0 < mvl; l0 += 8) {
            __m512i r[8];
            for (int t = 0; t < 8; t++) r[t] = GFNI_LD8(M + (size_t) (c + t) * mvl + l0);
            gfni_transpose_8x8_epi64(r);
            for (int l = 0; l < 8 && l0 + l < mvl; l++) GFNI_ST8(Mt + (size_t) (l0 + l) * GFNI_MT_S + c, r[l]);
        }
    }
    for (; c + 4 <= o; c += 4) {
        for (int l0 = 0; l0 < mvl; l0 += 4) {
            __m256i r[4];
            for (int t = 0; t < 4; t++) r[t] = GFNI_LD4(M + (size_t) (c + t) * mvl + l0);
            gfni_transpose_4x4_epi64(r);
            for (int l = 0; l < 4 && l0 + l < mvl; l++) GFNI_ST4(Mt + (size_t) (l0 + l) * GFNI_MT_S + c, r[l]);
        }
    }
    for (; c < o; c++) {   /* at most three columns */
        for (int l = 0; l < mvl; l++) Mt[(size_t) l * GFNI_MT_S + c] = M[(size_t) c * mvl + l];
    }
    for (int u = 0; u < nuse; u++) {
        const __m128i shl = _mm_cvtsi32_si128(bits[u]), shr = _mm_cvtsi32_si128(64 - bits[u]);
        int cc = 0;
        for (; cc + 8 <= o; cc += 8) GFNI_A_USE(__m512i, GFNI_LD8, GFNI_ST8, _mm512_xor_si512, _mm512_sll_epi64, _mm512_srl_epi64, GFNI_TERN8);
        for (; cc + 4 <= o; cc += 4) GFNI_A_USE(__m256i, GFNI_LD4, GFNI_ST4, _mm256_xor_si256, _mm256_sll_epi64, _mm256_srl_epi64, GFNI_TERN4);
        for (; cc + 2 <= o; cc += 2) GFNI_A_USE(__m128i, GFNI_LD2, GFNI_ST2, _mm_xor_si128, _mm_sll_epi64, _mm_srl_epi64, GFNI_TERN2);
        for (; cc < o; cc++) {
            uint64_t *d = A + (size_t) row0[u] * A_width + col[u] + cc;
            uint64_t prev = Mt[cc];
            d[0] ^= prev << bits[u];
            for (int l = 1; l < mvl; l++) {
                const uint64_t cur = Mt[(size_t) l * GFNI_MT_S + cc];
                d[(size_t) l * A_width] ^= bits[u] ? (cur << bits[u]) ^ (prev >> (64 - bits[u])) : cur;
                prev = cur;
            }
            if (bits[u] > 0) d[(size_t) mvl * A_width] ^= prev >> (64 - bits[u]);
        }
    }
    mayo_secure_clear(Mt, sizeof(Mt));
}
#endif
#endif
