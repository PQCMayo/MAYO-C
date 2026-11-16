// SPDX-License-Identifier: Apache-2.0

#ifndef SHUFFLE_ARITHMETIC_H
#define SHUFFLE_ARITHMETIC_H

#include <stdint.h>
#include <mayo.h>
#include <immintrin.h>
#include <string.h>
#include <stdalign.h>
#include <arithmetic_common.h>
#include <m_vec_arithmetic.h>

// P1/P2 are read straight from the packed keystream layout (m/2 bytes per m-vec),
// skipping the unpack spread. PINV maps an m-vec index to its packed byte address --
// the identity when m % 16 == 0. Index by number, not limb offset: no per-load
// divide, which these dual loads can't absorb.
#define MAYO_PACKED_P1P2 1
#if (M_MAX % 16) == 0
#define PINV(base, vec) ((const uint64_t *)(base) + (size_t)(vec) * M_VEC_LIMBS_MAX)
#else
#define PINV(base, vec) ((const uint64_t *)((const uint8_t *)(base) + (size_t)(vec) * (size_t)(M_MAX / 2)))
#endif


// multab storage, so callers can size buffers without knowing the engine
/* a 16-byte vpshufb table: low nibble = product for column 2q, high = 2q+1. Kernels
   broadcast it to both lanes; two adjacent tables load as [table q | table q+1]. */
typedef __m128i mayo_multab_t;
#define MAYO_V_MULTABS_N (((K_MAX + 1) / 2) * V_MAX)
#define MAYO_O_MULTABS_N ((O_AVX_ROUND_UP_ / 2) * V_MAX)
// select the multab sign orchestration in arithmetic.h
#define MAYO_HAVE_MULTAB_SIGN 1

/* Fused sign: expand_sk stores L' = P1*O + P2 instead of P2, and each attempt gets
   V*P1^t*O from the P1*V^t it computes anyway, skipping the P1^t*V^t sweep. L' is
   secret, so the caller clears it. Building L' costs about o columns of work, the
   sweep it saves about k: not for MAYO_2 (k = 5, o = 13). */
#if K_MAX >= O_MAX - 1
#define MAYO_SIGN_HALF_L 1
#endif

/* Products of a scalar matrix with a list of m-vecs. The m-vec is split into pieces
   that are each swept separately, so a piece's accumulators stay in registers:
   - groups of 4 limbs (one register). If the limb count is not a multiple of 4, the
     last group is loaded at limb M_VEC_LIMBS_MAX - 4 and its nibble mask zeroes the
     limbs an earlier group already covered;
   - if the limb count is 4n + 1 (MAYO_1, MAYO_5), the last limb is a piece of its own
     instead: its low and high nibbles fill one 128-bit lane, and the two lanes take
     two consecutive column pairs, so one shuffle covers four output columns. */
#define AVX_TAIL (M_VEC_LIMBS_MAX % 4 == 1 && M_VEC_LIMBS_MAX > 1)
#define AVX_GROUPS (AVX_TAIL ? M_VEC_LIMBS_MAX / 4 : (M_VEC_LIMBS_MAX + 3) / 4)
#define AVX_PIECES (AVX_GROUPS + AVX_TAIL)

static inline size_t avx_group_off(size_t g) {
    return 4 * g + 4 <= M_VEC_LIMBS_MAX ? 4 * g : M_VEC_LIMBS_MAX - 4;
}

static inline __m256i avx_group_mask(size_t g) {
    const int64_t lo = 0x0f0f0f0f0f0f0f0f;
    const size_t skip = 4 * g - avx_group_off(g); // limbs counted by group g - 1
    return _mm256_set_epi64x(lo, skip < 3 ? lo : 0, skip < 2 ? lo : 0, skip < 1 ? lo : 0);
}

// open piece g (a group, or the tail limb if g == AVX_GROUPS) and its accumulators
#define TILE_INIT(g, Len)                                                                \
    const int tail = (g) == AVX_GROUPS;                                                  \
    const size_t goff = avx_group_off(g);                                                \
    const __m256i gmask = avx_group_mask(g);                                             \
    __m256i temp[((Len) + 1) / 2 * 2] = {0}

/* temp += vec[piece] * multabs[0 .. Len). A group puts columns 2q, 2q+1 in temp[2q],
   temp[2q+1] (low and high nibbles of the m-vec bytes). The tail puts the tail limb's
   low and high nibbles in the two qwords of a lane and columns 4p .. 4p+3 in temp[p],
   reading tables q and q+1 with one 32-byte load. The last odd pair broadcasts table q
   instead, so that no load runs past the tables of the last source. */
#define TILE_MUL_ACC(vec, multabs, Len)                                                  \
    do {                                                                                 \
        if (!tail) {                                                                     \
            __m256i in_odd = _mm256_loadu_si256((const __m256i *) ((vec) + goff));       \
            __m256i in_even = _mm256_srli_epi16(in_odd, 4) & gmask;                      \
            in_odd &= gmask;                                                             \
            for (size_t k = 0; k < (Len); k += 2) {                                      \
                const __m256i tab = _mm256_broadcastsi128_si256((multabs)[k / 2]);      \
                temp[k]     ^= _mm256_shuffle_epi8(tab, in_odd);                         \
                temp[k + 1] ^= _mm256_shuffle_epi8(tab, in_even);                        \
            }                                                                            \
        } else {                                                                         \
            const __m256i x = _mm256_broadcastq_epi64(                                   \
                _mm_loadu_si64((const void *) ((vec) + M_VEC_LIMBS_MAX - 1)));           \
            const __m256i in = _mm256_srlv_epi64(x, _mm256_set_epi64x(4, 0, 4, 0))       \
                               & _mm256_set1_epi8(0x0f);                                 \
            for (size_t k = 0; k < (Len); k += 4) {                                      \
                const __m256i tab = k / 2 + 1 < ((Len) + 1) / 2                          \
                    ? _mm256_loadu_si256((const __m256i *) ((multabs) + k / 2))          \
                    : _mm256_broadcastsi128_si256((multabs)[k / 2]);                     \
                temp[k / 4] ^= _mm256_shuffle_epi8(tab, in);                             \
            }                                                                            \
        }                                                                                \
    } while (0)

/* acc[k * stride + piece] += column k of temp, k < Len. Per column pair the low-nibble
   products a and high-nibble products b of the source bytes become the two columns
   with t = (b ^ a >> 4) & 0x0f: a ^ t << 4 and b ^ t. The tail does the same per lane
   with a, b the two qwords, and xors the four resulting limbs into the last limb. */
#define TILE_ACC(Len, acc, stride)                                                       \
    do {                                                                                 \
        const __m256i lo = _mm256_set1_epi8(0x0f);                                       \
        if (!tail) {                                                                     \
            for (size_t k = 0; k < (Len); k += 2) {                                      \
                __m256i *a0 = (__m256i *) ((acc) + k * (stride) + goff);                 \
                __m256i t = (temp[k + 1] ^ _mm256_srli_epi16(temp[k], 4)) & lo;          \
                _mm256_storeu_si256(a0, _mm256_loadu_si256(a0) ^ temp[k] ^ _mm256_slli_epi16(t, 4)); \
                if (k + 1 < (Len)) {                                                     \
                    __m256i *a1 = (__m256i *) ((acc) + (k + 1) * (stride) + goff);       \
                    _mm256_storeu_si256(a1, _mm256_loadu_si256(a1) ^ temp[k + 1] ^ t);   \
                }                                                                        \
            }                                                                            \
        } else {                                                                         \
            for (size_t k = 0; k < (Len); k += 4) {                                      \
                const __m256i b = _mm256_unpackhi_epi64(temp[k / 4], temp[k / 4]);       \
                const __m256i t = (b ^ _mm256_srli_epi16(temp[k / 4], 4)) & lo;          \
                alignas(32) uint64_t cols[2][4];                                         \
                _mm256_store_si256((__m256i *) cols[0], temp[k / 4] ^ _mm256_slli_epi16(t, 4)); \
                _mm256_store_si256((__m256i *) cols[1], b ^ t);                          \
                for (size_t j = 0; j < 4 && k + j < (Len); j++) {                        \
                    (acc)[(k + j) * (stride) + M_VEC_LIMBS_MAX - 1] ^= cols[j % 2][2 * (j / 2)]; \
                }                                                                        \
            }                                                                            \
        }                                                                                \
    } while (0)

#define MVL M_VEC_LIMBS_MAX
#define O_TABS(r) (O_multabs + O_AVX_ROUND_UP_ / 2 * (r))
#define K_TABS(tabs, r) ((tabs) + K_OVER_2 * (r))

// P1 * O -> P2: P1 is v x v, O is v x o
static
inline void P1_times_O(const uint64_t *P1, mayo_multab_t *O_multabs, uint64_t *acc){
    size_t row_start = 0;
    for (size_t r = 0; r < V_MAX; r++) {
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, O_MAX);
            for (size_t c = r; c < V_MAX; c++) {
                TILE_MUL_ACC(PINV(P1, row_start + c - r), O_TABS(c), O_MAX);
            }
            TILE_ACC(O_MAX, acc + r*O_MAX*MVL, MVL);
        }
        row_start += V_MAX - r;
    }
}

static
inline void Ot_times_P1O_P2(const uint64_t *P1O_P2, mayo_multab_t *O_multabs, uint64_t *acc){
    for (size_t c = 0; c < O_MAX; c++) {
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, O_MAX);
            for (size_t r = 0; r < V_MAX; r++) {
                TILE_MUL_ACC(P1O_P2 + (r*O_MAX + c)*MVL, O_TABS(r), O_MAX);
            }
            TILE_ACC(O_MAX, acc + c*MVL, O_MAX*MVL);
        }
    }
}

static
inline void P1P1t_times_O(const mayo_params_t* p, const uint64_t *P1, const unsigned char *O, uint64_t *acc){
    (void) p;
    mayo_multab_t O_multabs[MAYO_O_MULTABS_N];
    mayo_O_multabs(O, O_multabs);

    size_t row_start = 0;
    for (size_t r = 0; r < V_MAX; r++) {
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, O_MAX);
            size_t pos = r; // P1[c][r], column r above the diagonal
            for (size_t c = 0; c < r; c++) {
                TILE_MUL_ACC(PINV(P1, pos), O_TABS(c), O_MAX);
                pos += V_MAX - c - 1;
            }
            for (size_t c = r + 1; c < V_MAX; c++) { // P1[r][c], row r right of the diagonal
                TILE_MUL_ACC(PINV(P1, row_start + c - r), O_TABS(c), O_MAX);
            }
            TILE_ACC(O_MAX, acc + r*O_MAX*MVL, MVL);
        }
        row_start += V_MAX - r;
    }
    mayo_secure_clear(O_multabs, sizeof(O_multabs));
}

/* acc[c] += sum_r O[r][c] * lambda[r]: the one column of U^t * O that the
   materialized-L signing path has to supply itself. Ut_times_O below with its
   column loop removed and lambda read in place of a column of U. */
#define MAYO_HAVE_LAMBDA_TIMES_O 1
static inline void m_lambda_times_O(const uint64_t *lambda, const mayo_multab_t *O_multabs,
                                    uint64_t *acc) {
    for (size_t g = 0; g < AVX_PIECES; g++) {
        TILE_INIT(g, O_MAX);
        for (size_t r = 0; r < V_MAX; r++) {
            TILE_MUL_ACC(lambda + r*MVL, O_TABS(r), O_MAX);
        }
        TILE_ACC(O_MAX, acc, MVL);
    }
}

// acc += U^t * O, where U is a v-by-k matrix of m-vectors.
static
inline void Ut_times_O(const uint64_t *U, const mayo_multab_t *O_multabs, uint64_t *acc){
    for (size_t c = 0; c < K_MAX; c++) {
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, O_MAX);
            for (size_t r = 0; r < V_MAX; r++) {
                TILE_MUL_ACC(U + (r*K_MAX + c)*MVL, O_TABS(r), O_MAX);
            }
            TILE_ACC(O_MAX, acc + c*O_MAX*MVL, MVL);
        }
    }
}

static
inline void Vt_times_L(const uint64_t *L, const mayo_multab_t *V_multabs, uint64_t *acc){
    for (size_t c = 0; c < O_MAX; c++) {
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, K_MAX);
            for (size_t r = 0; r < V_MAX; r++) {
                TILE_MUL_ACC(L + (r*O_MAX + c)*MVL, K_TABS(V_multabs, r), K_MAX);
            }
            TILE_ACC(K_MAX, acc + c*MVL, O_MAX*MVL);
        }
    }
}

static
inline void Vt_times_Pv(const uint64_t *Pv, const mayo_multab_t *V_multabs, uint64_t *acc){
    for (size_t c = 0; c < K_MAX; c++) {
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, K_MAX);
            for (size_t r = 0; r < V_MAX; r++) {
                TILE_MUL_ACC(Pv + (r*K_MAX + c)*MVL, K_TABS(V_multabs, r), K_MAX);
            }
            TILE_ACC(K_MAX, acc + c*MVL, K_MAX*MVL);
        }
    }
}

static
inline void P1_times_Vt(const uint64_t *P1, mayo_multab_t *V_multabs, uint64_t *acc){
    size_t row_start = 0;
    for (size_t r = 0; r < V_MAX; r++) {
        memset(acc + r*K_MAX*MVL, 0, K_MAX*MVL*sizeof(uint64_t)); // acc = P1 * V^t
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, K_MAX);
            for (size_t c = r; c < V_MAX; c++) {
                TILE_MUL_ACC(PINV(P1, row_start + c - r), K_TABS(V_multabs, c), K_MAX);
            }
            TILE_ACC(K_MAX, acc + r*K_MAX*MVL, MVL);
        }
        row_start += V_MAX - r;
    }
}

// acc += P1^t * V^t, where P1 is stored as an upper-triangular matrix.
static
inline void P1t_times_Vt(const uint64_t *P1, mayo_multab_t *V_multabs, uint64_t *acc){
    size_t row_start = 0;
    for (size_t r = 0; r < V_MAX; r++) {
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, K_MAX);
            size_t pos = r; // P1[c][r] for c <= r
            for (size_t c = 0; c < r; c++) {
                TILE_MUL_ACC(PINV(P1, pos), K_TABS(V_multabs, c), K_MAX);
                pos += V_MAX - c - 1;
            }
            TILE_MUL_ACC(PINV(P1, row_start), K_TABS(V_multabs, r), K_MAX);
            TILE_ACC(K_MAX, acc + r*K_MAX*MVL, MVL);
        }
        row_start += V_MAX - r;
    }
}

// S1^t * PS1: the same product as V * Pv
static
inline void S1t_times_PS1(const uint64_t *_PS1, mayo_multab_t *S1_multabs, uint64_t *_acc){
    Vt_times_Pv(_PS1, S1_multabs, _acc);
}

static
inline void S2t_times_PS2(const uint64_t *PS2, mayo_multab_t *S2_multabs, uint64_t *acc){
    for (size_t c = 0; c < K_MAX; c++) {
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, K_MAX);
            for (size_t r = 0; r < O_MAX; r++) {
                TILE_MUL_ACC(PS2 + (r*K_MAX + c)*MVL, K_TABS(S2_multabs, r), K_MAX);
            }
            TILE_ACC(K_MAX, acc + c*MVL, K_MAX*MVL);
        }
    }
}

// PS rows 0..v-1: P1 * S1 + P2 * S2
static
inline void P1_times_S1_plus_P2_times_S2(const uint64_t *P1, const uint64_t *P2, mayo_multab_t *S1_multabs, mayo_multab_t *S2_multabs, uint64_t *acc){
    size_t row_start = 0;
    for (size_t r = 0; r < V_MAX; r++) {
        memset(acc + r*K_MAX*MVL, 0, K_MAX*MVL*sizeof(uint64_t)); // acc = P1*S1 + P2*S2
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, K_MAX);
            for (size_t c = r; c < V_MAX; c++) { // P1 * S1
                TILE_MUL_ACC(PINV(P1, row_start + c - r), K_TABS(S1_multabs, c), K_MAX);
            }
            for (size_t c = 0; c < O_MAX; c++) { // P2 * S2
                TILE_MUL_ACC(PINV(P2, r*O_MAX + c), K_TABS(S2_multabs, c), K_MAX);
            }
            TILE_ACC(K_MAX, acc + r*K_MAX*MVL, MVL);
        }
        row_start += V_MAX - r;
    }
}

// P3*S2 -> P3: o x o, S2: o x k // P3 upper triangular
static
inline void P3_times_S2(const uint64_t *P3, mayo_multab_t *S2_multabs, uint64_t *acc){
    size_t row_start = 0;
    for (size_t r = 0; r < O_MAX; r++) {
        memset(acc + r*K_MAX*MVL, 0, K_MAX*MVL*sizeof(uint64_t)); // acc = P3 * S2
        for (size_t g = 0; g < AVX_PIECES; g++) {
            TILE_INIT(g, K_MAX);
            for (size_t c = r; c < O_MAX; c++) {
                TILE_MUL_ACC(P3 + (row_start + c - r)*MVL, K_TABS(S2_multabs, c), K_MAX);
            }
            TILE_ACC(K_MAX, acc + r*K_MAX*MVL, MVL);
        }
        row_start += O_MAX - r;
    }
}


static inline 
void compute_P3(const mayo_params_t* p, const uint64_t* P1, uint64_t *P2, const unsigned char *O, uint64_t *P3){
    (void) p;
    mayo_multab_t O_multabs[MAYO_O_MULTABS_N];
    mayo_O_multabs(O, O_multabs);
    P1_times_O(P1, O_multabs, P2);
    Ot_times_P1O_P2(P2, O_multabs, P3);
    mayo_secure_clear(O_multabs, sizeof(O_multabs));
}

// P2 <- L' = P1*O + P2, the first half of compute_P3
static inline
void compute_L_half(const uint64_t *P1, uint64_t *P2, const unsigned char *O){
    mayo_multab_t O_multabs[MAYO_O_MULTABS_N];
    mayo_O_multabs(O, O_multabs);
    P1_times_O(P1, O_multabs, P2);
    mayo_secure_clear(O_multabs, sizeof(O_multabs));
}

// compute P * S^t = [ P1  P2 ] * [S1] = [P1*S1 + P2*S2]
//                   [  0  P3 ]   [S2]   [        P3*S2]
// compute S * PS  = [ S1 S2 ] * [ P1*S1 + P2*S2 = P1 ] = [ S1*P1 + S2*P2 ]
//                               [         P3*S2 = P2 ]
static inline void m_calculate_PS_SPS(const mayo_params_t *p, const uint64_t *P1, const uint64_t *P2, const uint64_t *P3, const unsigned char *S,
                              const uint64_t *lambda, uint64_t *SPS) {
    (void) p;
    const int o = PARAM_NAME(o);
    const int v = PARAM_NAME(v);
    const int k = PARAM_NAME(k);
    const int n = o + v;
    unsigned char S1[V_MAX*K_MAX]; // == N-O, K
    unsigned char S2[O_MAX*K_MAX]; // == O, K
    for (int r = 0; r < k; r++) {
        memcpy(S1 + r*v, S + r*n, v);
        memcpy(S2 + r*o, S + r*n + v, o);
    }

    alignas(64) uint64_t PS[N_MAX * K_MAX * M_VEC_LIMBS_MAX]; // each row zeroed by its kernel

    mayo_multab_t S1_multabs[(K_MAX+1)/2*V_MAX];    
    mayo_multab_t S2_multabs[(K_MAX+1)/2*O_MAX];
    mayo_S_multabs(S1, V_MAX, S1_multabs);
    mayo_S_multabs(S2, O_MAX, S2_multabs);

    P1_times_S1_plus_P2_times_S2(P1, P2, S1_multabs, S2_multabs, PS);
    P3_times_S2(P3, S2_multabs, PS + V_MAX*K_MAX*M_VEC_LIMBS_MAX); // upper triangular

    /* lambda occupies column 0 of PS. The S^t * PS pass below applies S[i][r] to it,
       so SPS[i][0] picks up Lambda(s_i) for n m-vector XORs and no multiplications. */
    for (int r = 0; r < v; r++) {
        uint64_t *dst = PS + (size_t) r * K_MAX * M_VEC_LIMBS_MAX;
        for (int l = 0; l < M_VEC_LIMBS_MAX; l++) {
            dst[l] ^= lambda[(size_t) r * M_VEC_LIMBS_MAX + l];
        }
    }
    for (int r = 0; r < o; r++) {
        uint64_t *dst = PS + V_MAX * K_MAX * M_VEC_LIMBS_MAX + (size_t) r * K_MAX * M_VEC_LIMBS_MAX;
        for (int l = 0; l < M_VEC_LIMBS_MAX; l++) {
            dst[l] ^= lambda[(size_t) (v + r) * M_VEC_LIMBS_MAX + l];
        }
    }

    // S^T * PS = S1^t*PS1 + S2^t*PS2
    S1t_times_PS1(PS, S1_multabs, SPS);
    S2t_times_PS2(PS + V_MAX*K_MAX*M_VEC_LIMBS_MAX, S2_multabs, SPS);
}

#undef K_OVER_2

/* 4x4 and 2x2 transposes of 64-bit lanes */
static inline void avx2_transpose_4x4_epi64(__m256i *r) {
    __m256i t0 = _mm256_unpacklo_epi64(r[0], r[1]);
    __m256i t1 = _mm256_unpackhi_epi64(r[0], r[1]);
    __m256i t2 = _mm256_unpacklo_epi64(r[2], r[3]);
    __m256i t3 = _mm256_unpackhi_epi64(r[2], r[3]);
    r[0] = _mm256_permute2x128_si256(t0, t2, 0x20);
    r[1] = _mm256_permute2x128_si256(t1, t3, 0x20);
    r[2] = _mm256_permute2x128_si256(t0, t2, 0x31);
    r[3] = _mm256_permute2x128_si256(t1, t3, 0x31);
}

/* see the fallback in arithmetic.h. Whole blocks of columns are transposed in registers
   and xored into A rows; only whole ones, since a masked store costs a full store either
   way. The loads run past the last column into Mtmp's MAYO_MVEC_SLACK tail. */
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
            avx2_transpose_4x4_epi64(r);
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
            const __m128i rr[2] = { _mm_unpacklo_epi64(a0, a1), _mm_unpackhi_epi64(a0, a1) };
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


/* dst[i] ^= each nibble of src[i] times coeff, one shuffle per nibble half */
#define MAYO_HAVE_GF16_MULC_XOR 1
#define MAYO_GF16_MULC_NEEDS_TABLE 1
static inline void m_gf16_mulc_xor(uint64_t *dst, const uint64_t *src, int n, unsigned char coeff, const unsigned char *t16) {
    (void) coeff;
    const __m256i tab = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *) t16));
    const __m256i lown = _mm256_set1_epi8(0x0f);
    const int batched = n & ~3;
    for (int i = 0; i < batched; i += 4) {
        __m256i v = _mm256_loadu_si256((const __m256i *) (src + i));
        __m256i lo = _mm256_and_si256(v, lown);
        __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), lown);
        __m256i r = _mm256_xor_si256(_mm256_shuffle_epi8(tab, lo),
                                     _mm256_slli_epi16(_mm256_shuffle_epi8(tab, hi), 4));
        _mm256_storeu_si256((__m256i *) (dst + i),
            _mm256_xor_si256(_mm256_loadu_si256((const __m256i *) (dst + i)), r));
    }
    for (int i = batched; i < n; i++) {
        const uint64_t low = 0x1111111111111111ULL, v = src[i];
        dst[i] ^= ( v       & low) * t16[1] ^ ((v >> 1) & low) * t16[2]
                ^ ((v >> 2) & low) * t16[4] ^ ((v >> 3) & low) * t16[8];
    }
}

#endif
