// SPDX-License-Identifier: Apache-2.0 and MIT and Public Domain

#ifdef ENABLE_AESNI

#include <mem.h>
#include <stdint.h>
#include <string.h>
#include <tmmintrin.h>
#include <wmmintrin.h>
#if defined(__VAES__)
#include <immintrin.h>
#endif

// Adapted from liboqs/src/common/aes which in turn takes it from:
// crypto_core/aes128ncrypt/dolbeau/aesenc-int
// (https://bench.cr.yp.to/supercop.html)
static inline void aes128ni_setkey_encrypt(const unsigned char *key,
        __m128i rkeys[11]) {
    __m128i key0 = _mm_loadu_si128((const __m128i *)(key + 0));
    __m128i temp0, temp1, temp4;
    int idx = 0;

    temp0 = key0;

#define BLOCK1(IMM)                                                            \
  temp1 = _mm_aeskeygenassist_si128(temp0, IMM);                               \
  rkeys[idx++] = temp0;                                                        \
  temp4 = _mm_slli_si128(temp0, 4);                                            \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp4 = _mm_slli_si128(temp0, 8);                                            \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp1 = _mm_shuffle_epi32(temp1, 0xff);                                      \
  temp0 = _mm_xor_si128(temp0, temp1)

    BLOCK1(0x01);
    BLOCK1(0x02);
    BLOCK1(0x04);
    BLOCK1(0x08);
    BLOCK1(0x10);
    BLOCK1(0x20);
    BLOCK1(0x40);
    BLOCK1(0x80);
    BLOCK1(0x1b);
    BLOCK1(0x36);
    rkeys[idx++] = temp0;
}

void oqs_aes128_load_schedule_ni(const uint8_t *key, void **_schedule) {
    *_schedule = malloc(11 * sizeof(__m128i));
    // assert(*_schedule != NULL);
    __m128i *schedule = (__m128i *)*_schedule;
    aes128ni_setkey_encrypt(key, schedule);
}

void oqs_aes128_free_schedule_ni(void *schedule) {
    if (schedule != NULL) {
        mayo_secure_free(schedule, 11 * sizeof(__m128i));
    }
}

// Single encryption
static inline void aes128ni_encrypt(const __m128i rkeys[11], __m128i nv,
                                    unsigned char *out) {
    __m128i temp = _mm_xor_si128(nv, rkeys[0]);
    temp = _mm_aesenc_si128(temp, rkeys[1]);
    temp = _mm_aesenc_si128(temp, rkeys[2]);
    temp = _mm_aesenc_si128(temp, rkeys[3]);
    temp = _mm_aesenc_si128(temp, rkeys[4]);
    temp = _mm_aesenc_si128(temp, rkeys[5]);
    temp = _mm_aesenc_si128(temp, rkeys[6]);
    temp = _mm_aesenc_si128(temp, rkeys[7]);
    temp = _mm_aesenc_si128(temp, rkeys[8]);
    temp = _mm_aesenc_si128(temp, rkeys[9]);
    temp = _mm_aesenclast_si128(temp, rkeys[10]);
    _mm_storeu_si128((__m128i *)(out), temp);
}

// 4x interleaved encryption
static inline void aes128ni_encrypt_x4(const __m128i rkeys[11], __m128i n0,
                                       __m128i n1, __m128i n2, __m128i n3,
                                       unsigned char *out) {
    __m128i temp0 = _mm_xor_si128(n0, rkeys[0]);
    __m128i temp1 = _mm_xor_si128(n1, rkeys[0]);
    __m128i temp2 = _mm_xor_si128(n2, rkeys[0]);
    __m128i temp3 = _mm_xor_si128(n3, rkeys[0]);

#define AESNENCX4(IDX)                                                         \
  temp0 = _mm_aesenc_si128(temp0, rkeys[IDX]);                                 \
  temp1 = _mm_aesenc_si128(temp1, rkeys[IDX]);                                 \
  temp2 = _mm_aesenc_si128(temp2, rkeys[IDX]);                                 \
  temp3 = _mm_aesenc_si128(temp3, rkeys[IDX])

    AESNENCX4(1);
    AESNENCX4(2);
    AESNENCX4(3);
    AESNENCX4(4);
    AESNENCX4(5);
    AESNENCX4(6);
    AESNENCX4(7);
    AESNENCX4(8);
    AESNENCX4(9);

    temp0 = _mm_aesenclast_si128(temp0, rkeys[10]);
    temp1 = _mm_aesenclast_si128(temp1, rkeys[10]);
    temp2 = _mm_aesenclast_si128(temp2, rkeys[10]);
    temp3 = _mm_aesenclast_si128(temp3, rkeys[10]);

    _mm_storeu_si128((__m128i *)(out + 0), temp0);
    _mm_storeu_si128((__m128i *)(out + 16), temp1);
    _mm_storeu_si128((__m128i *)(out + 32), temp2);
    _mm_storeu_si128((__m128i *)(out + 48), temp3);
}

// VAES bulk CTR: 16 blocks/iter (VAES-512) or 8 (VAES-256); the counter is
// byte-reversed to big-endian just before AES, matching the scalar tail.
// VAES-512 is gated on the gfni backend: 512-bit ops throttle nearby 256-bit
// vpshufb throughput, so it would only slow the avx2 arithmetic, with no AES
// gain over VAES-256.
#if defined(__VAES__) && defined(__AVX512F__) && defined(MAYO_GFNI)
static inline void aes_ctr_vaes_bulk(const __m128i *rkeys, int nrounds, uint8_t **outp,
                                     size_t *lenp, __m128i *block, __m128i mask) {
    uint8_t *out = *outp;
    size_t out_len = *lenp;
    __m512i brk[11];
    for (int i = 0; i <= nrounds; i++) brk[i] = _mm512_broadcast_i32x4(rkeys[i]);
    const __m512i M = _mm512_broadcast_i32x4(mask);
    const __m512i INC = _mm512_set_epi64(16, 0, 16, 0, 16, 0, 16, 0);
    const __m512i P4 = _mm512_set_epi64(4, 0, 4, 0, 4, 0, 4, 0);
    __m512i c[4];
    c[0] = _mm512_set_epi64(3, 0, 2, 0, 1, 0, 0, 0);
    for (int j = 1; j < 4; j++) c[j] = _mm512_add_epi64(c[j - 1], P4);
    uint64_t count = 0;
    while (out_len >= 256) {
        __m512i t[4];
        for (int j = 0; j < 4; j++) t[j] = _mm512_xor_si512(_mm512_shuffle_epi8(c[j], M), brk[0]);
        for (int r = 1; r < nrounds; r++)
            for (int j = 0; j < 4; j++) t[j] = _mm512_aesenc_epi128(t[j], brk[r]);
        for (int j = 0; j < 4; j++) {
            t[j] = _mm512_aesenclast_epi128(t[j], brk[nrounds]);
            _mm512_storeu_si512((void *)(out + 64 * j), t[j]);
            c[j] = _mm512_add_epi64(c[j], INC);
        }
        out += 256; out_len -= 256; count += 16;
    }
    *outp = out;
    *lenp = out_len;
    *block = _mm_shuffle_epi8(_mm_set_epi64x((long long)count, 0), mask);
}
#define AES_HAVE_VAES 1
#elif defined(__VAES__)
static inline void aes_ctr_vaes_bulk(const __m128i *rkeys, int nrounds, uint8_t **outp,
                                     size_t *lenp, __m128i *block, __m128i mask) {
    uint8_t *out = *outp;
    size_t out_len = *lenp;
    __m256i brk[11];
    for (int i = 0; i <= nrounds; i++) brk[i] = _mm256_broadcastsi128_si256(rkeys[i]);
    const __m256i M = _mm256_broadcastsi128_si256(mask);
    const __m256i INC = _mm256_set_epi64x(8, 0, 8, 0);
    const __m256i P2 = _mm256_set_epi64x(2, 0, 2, 0);
    __m256i c[4];
    c[0] = _mm256_set_epi64x(1, 0, 0, 0);
    for (int j = 1; j < 4; j++) c[j] = _mm256_add_epi64(c[j - 1], P2);
    uint64_t count = 0;
    while (out_len >= 128) {
        __m256i t[4];
        for (int j = 0; j < 4; j++) t[j] = _mm256_xor_si256(_mm256_shuffle_epi8(c[j], M), brk[0]);
        for (int r = 1; r < nrounds; r++)
            for (int j = 0; j < 4; j++) t[j] = _mm256_aesenc_epi128(t[j], brk[r]);
        for (int j = 0; j < 4; j++) {
            t[j] = _mm256_aesenclast_epi128(t[j], brk[nrounds]);
            _mm256_storeu_si256((__m256i *)(out + 32 * j), t[j]);
            c[j] = _mm256_add_epi64(c[j], INC);
        }
        out += 128; out_len -= 128; count += 8;
    }
    *outp = out;
    *lenp = out_len;
    *block = _mm_shuffle_epi8(_mm_set_epi64x((long long)count, 0), mask);
}
#define AES_HAVE_VAES 1
#endif

// Not for general use: IV = 0, nonce = 0
static void oqs_aes128_ctr_enc_sch_ni(const void *schedule, uint8_t *out,
                                      size_t out_len) {
    __m128i mask =
        _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
    __m128i block = _mm_set_epi64x(0, 0);
    // block = _mm_xor_si128(block, block); // set to zero
#ifdef AES_HAVE_VAES
    aes_ctr_vaes_bulk((const __m128i *)schedule, 10, &out, &out_len, &block, mask);
#endif

    while (out_len >= 64) {
        __m128i nv0 = block;
        __m128i nv1 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)),
                          mask);
        __m128i nv2 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(2, 0)),
                          mask);
        __m128i nv3 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(3, 0)),
                          mask);
        aes128ni_encrypt_x4(schedule, nv0, nv1, nv2, nv3, out);
        block = _mm_shuffle_epi8(
                    _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(4, 0)),
                    mask);
        out += 64;
        out_len -= 64;
    }
    while (out_len >= 16) {
        aes128ni_encrypt(schedule, block, out);
        out += 16;
        out_len -= 16;
        block = _mm_shuffle_epi8(
                    _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)),
                    mask);
    }
    if (out_len > 0) {
        uint8_t tmp[16];
        aes128ni_encrypt(schedule, block, tmp);
        memcpy(out, tmp, out_len);
    }
}

/* AES-128 ECB over nblocks whole blocks; matches AES_128_ECB_C exactly. */
int AES_128_ECB_NI(unsigned char *output, const unsigned char *input, size_t nblocks,
                   const unsigned char *key) {
    __m128i rkeys[11];
    aes128ni_setkey_encrypt(key, rkeys);
    size_t i = 0;
    for (; i + 4 <= nblocks; i += 4) {
        aes128ni_encrypt_x4(rkeys,
                            _mm_loadu_si128((const __m128i *)(input + 16 * (i + 0))),
                            _mm_loadu_si128((const __m128i *)(input + 16 * (i + 1))),
                            _mm_loadu_si128((const __m128i *)(input + 16 * (i + 2))),
                            _mm_loadu_si128((const __m128i *)(input + 16 * (i + 3))),
                            output + 16 * i);
    }
    for (; i < nblocks; i++) {
        aes128ni_encrypt(rkeys, _mm_loadu_si128((const __m128i *)(input + 16 * i)),
                         output + 16 * i);
    }
    mayo_secure_clear(rkeys, sizeof rkeys);
    return (int)(nblocks * 16);
}

int AES_128_CTR_NI(unsigned char *output, size_t outputByteLen,
                   const unsigned char *input, size_t inputByteLen) {
    void *schedule = NULL;
    oqs_aes128_load_schedule_ni(input, &schedule);
    oqs_aes128_ctr_enc_sch_ni(schedule, output, outputByteLen);
    oqs_aes128_free_schedule_ni(schedule);
    return (int)outputByteLen;
}

// 4-Round AES...

// From crypto_core/aes128ncrypt/dolbeau/aesenc-int
static inline void aes128r4ni_setkey_encrypt(const unsigned char *key,
        __m128i rkeys[5]) {
    __m128i key0 = _mm_loadu_si128((const __m128i *)(key + 0));
    __m128i temp0, temp1, temp4;
    int idx = 0;

    temp0 = key0;

    /* blockshift-based block by Cedric Bourrasset */
#define BLOCK1(IMM)                                                            \
  temp1 = _mm_aeskeygenassist_si128(temp0, IMM);                               \
  rkeys[idx++] = temp0;                                                        \
  temp4 = _mm_slli_si128(temp0, 4);                                            \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp4 = _mm_slli_si128(temp0, 8);                                            \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp1 = _mm_shuffle_epi32(temp1, 0xff);                                      \
  temp0 = _mm_xor_si128(temp0, temp1)

    BLOCK1(0x01);
    BLOCK1(0x02);
    BLOCK1(0x04);
    BLOCK1(0x08);
    rkeys[idx++] = temp0;
}

void oqs_aes128r4_load_schedule_ni(const uint8_t *key, void **_schedule) {
    *_schedule = malloc(5 * sizeof(__m128i));
    // assert(*_schedule != NULL);
    __m128i *schedule = (__m128i *)*_schedule;
    aes128r4ni_setkey_encrypt(key, schedule);
}

void oqs_aes128r4_free_schedule_ni(void *schedule) {
    if (schedule != NULL) {
        mayo_secure_free(schedule, 5 * sizeof(__m128i));
    }
}

// Single encryption
static inline void aes128r4ni_encrypt(const __m128i rkeys[5], __m128i nv,
                                      unsigned char *out) {
    __m128i temp = _mm_xor_si128(nv, rkeys[0]);
    temp = _mm_aesenc_si128(temp, rkeys[1]);
    temp = _mm_aesenc_si128(temp, rkeys[2]);
    temp = _mm_aesenc_si128(temp, rkeys[3]);
    temp = _mm_aesenclast_si128(temp, rkeys[4]);
    _mm_storeu_si128((__m128i *)(out), temp);
}

// 4x interleaved encryption
static inline void aes128r4ni_encrypt_x4(const __m128i rkeys[5], __m128i n0,
        __m128i n1, __m128i n2, __m128i n3,
        unsigned char *out) {
    __m128i temp0 = _mm_xor_si128(n0, rkeys[0]);
    __m128i temp1 = _mm_xor_si128(n1, rkeys[0]);
    __m128i temp2 = _mm_xor_si128(n2, rkeys[0]);
    __m128i temp3 = _mm_xor_si128(n3, rkeys[0]);

#define AESNENCX4(IDX)                                                         \
  temp0 = _mm_aesenc_si128(temp0, rkeys[IDX]);                                 \
  temp1 = _mm_aesenc_si128(temp1, rkeys[IDX]);                                 \
  temp2 = _mm_aesenc_si128(temp2, rkeys[IDX]);                                 \
  temp3 = _mm_aesenc_si128(temp3, rkeys[IDX])

    AESNENCX4(1);
    AESNENCX4(2);
    AESNENCX4(3);

    temp0 = _mm_aesenclast_si128(temp0, rkeys[4]);
    temp1 = _mm_aesenclast_si128(temp1, rkeys[4]);
    temp2 = _mm_aesenclast_si128(temp2, rkeys[4]);
    temp3 = _mm_aesenclast_si128(temp3, rkeys[4]);

    _mm_storeu_si128((__m128i *)(out + 0), temp0);
    _mm_storeu_si128((__m128i *)(out + 16), temp1);
    _mm_storeu_si128((__m128i *)(out + 32), temp2);
    _mm_storeu_si128((__m128i *)(out + 48), temp3);
}

// Not for general use: IV = 0, nonce = 0
static void oqs_aes128r4_ctr_enc_sch_ni(const void *schedule, uint8_t *out,
                                        size_t out_len) {
    __m128i mask =
        _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
    __m128i block = _mm_set_epi64x(0, 0);
#ifdef AES_HAVE_VAES
    aes_ctr_vaes_bulk((const __m128i *)schedule, 4, &out, &out_len, &block, mask);
#endif

    while (out_len >= 64) {
        __m128i nv0 = block;
        __m128i nv1 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)),
                          mask);
        __m128i nv2 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(2, 0)),
                          mask);
        __m128i nv3 = _mm_shuffle_epi8(
                          _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(3, 0)),
                          mask);
        aes128r4ni_encrypt_x4(schedule, nv0, nv1, nv2, nv3, out);
        block = _mm_shuffle_epi8(
                    _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(4, 0)),
                    mask);
        out += 64;
        out_len -= 64;
    }
    while (out_len >= 16) {
        aes128r4ni_encrypt(schedule, block, out);
        out += 16;
        out_len -= 16;
        block = _mm_shuffle_epi8(
                    _mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)),
                    mask);
    }
    if (out_len > 0) {
        uint8_t tmp[16];
        aes128r4ni_encrypt(schedule, block, tmp);
        memcpy(out, tmp, out_len);
    }
}

int AES_128_CTR_4R_NI(unsigned char *output, size_t outputByteLen,
                      const unsigned char *input, size_t inputByteLen) {
    void *schedule = NULL;
    oqs_aes128r4_load_schedule_ni(input, &schedule);
    oqs_aes128r4_ctr_enc_sch_ni(schedule, output, outputByteLen);
    oqs_aes128r4_free_schedule_ni(schedule);
    return (int)outputByteLen;
}
#endif

