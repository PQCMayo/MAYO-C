// SPDX-License-Identifier: MIT and Public Domain

#ifdef ENABLE_AESNEON

// Code taken from https://github.com/ChristerKnorborg/post-quantum-signature-schemes/blob/main/src/genkat/aes_arm.c

/* ARMv8 AES Implementation adapted
  from liboqs/src/common/aes which
  in turn takes it from:
  crypto_core/aes128ncrypt/dolbeau/aesenc-int
 (https://bench.cr.yp.to/supercop.html) */

#include "mem.h"
#include <string.h>
#include <stdlib.h>

#if defined(__arm__) || defined(__aarch32__) || defined(__arm64__) || defined(__aarch64__) || defined(_M_ARM) || defined(_M_ARM64)
# if defined(__GNUC__)
#  include <stdint.h>
# endif
# if defined(__ARM_NEON) || defined(_MSC_VER)
#  include <arm_neon.h>
# endif
/* GCC and LLVM Clang, but not Apple Clang */
# if defined(__GNUC__) && !defined(__apple_build_version__)
#  if defined(__ARM_ACLE) || defined(__ARM_FEATURE_CRYPTO)
#   include <arm_acle.h>
#  endif
# endif
#endif  /* ARM Headers */

// aes s-box
static const uint8_t sbox[256] = {
  //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };




// subword algorithm used in the aes key scheduling.
uint32_t subword(uint32_t word) {
    return (uint32_t)sbox[(word >> 24) & 0xFF] << 24 |
           (uint32_t)sbox[(word >> 16) & 0xFF] << 16 |
           (uint32_t)sbox[(word >> 8) & 0xFF] << 8 |
           (uint32_t)sbox[word & 0xFF];
}

uint32x4_t aeskeygenassist(uint32x4_t a32, uint8_t rcon) {
    // Extract words X1 and X3
    uint32_t X1 = vgetq_lane_u32(a32, 1);
    uint32_t X3 = vgetq_lane_u32(a32, 3);

    // Apply SubWord (Manually implemented using AES S-box lookup table)
    uint32_t subX1 = subword(X1);  // Implement this function based on AES S-box
    uint32_t subX3 = subword(X3);  // Implement this function based on AES S-box

    // RotWord
    uint32_t rotX1 = (subX1 >> 8) | (subX1 << 24);
    uint32_t rotX3 = (subX3 >> 8) | (subX3 << 24);

    // Apply RCON
    rotX1 ^= (uint32_t)rcon;
    rotX3 ^= (uint32_t)rcon;

    // Assemble the final vector
    uint32x4_t result = a32;
    result = vsetq_lane_u32(subX1, result, 0);
    result = vsetq_lane_u32(rotX1, result, 1);
    result = vsetq_lane_u32(subX3, result, 2);
    result = vsetq_lane_u32(rotX3, result, 3);

    return result;
}

static void aes_setkey_encrypt(const unsigned char *key, uint8x16_t rkeys[]) {
    uint8x16_t key0 = vld1q_u8((const uint8_t *)(key));
    uint32x4_t temp0, temp1, temp4;
    int idx = 0;

    temp0 = vreinterpretq_u32_u8(key0);
    temp4 = vdupq_n_u32(0);

    #define BLOCK1(IMM)                                                     \
      temp1 = aeskeygenassist(temp0, IMM);                                  \
      rkeys[idx++] = vreinterpretq_u8_u32(temp0);                           \
      temp4 = vsetq_lane_u32(vgetq_lane_u32(temp0, 0), temp4, 1);           \
      temp4 = vsetq_lane_u32(vgetq_lane_u32(temp0, 1), temp4, 2);           \
      temp4 = vsetq_lane_u32(vgetq_lane_u32(temp0, 2), temp4, 3);           \
      temp0 = veorq_u32(temp0, temp4);                                      \
      temp0 = vreinterpretq_u32_u64((vsetq_lane_u64(((uint64_t) veor_u32(vget_high_u32(temp0), vget_low_u32(temp0))), vreinterpretq_u64_u32(temp0), 1)));      \
      temp1 = vdupq_n_u32(vgetq_lane_u32(temp1, 3));                        \
      temp0 = veorq_u32(temp0, temp1);                                      \

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
    rkeys[idx++] = vreinterpretq_u8_u32(temp0);
}

void arm_aes128_load_schedule(const uint8_t *key, void **_schedule) {
    *_schedule = malloc(11 * sizeof(uint8x16_t));
    // assert(*_schedule != NULL);
    uint8x16_t *schedule = (uint8x16_t *)*_schedule;
    aes_setkey_encrypt(key, schedule);
}


// AES encryption using NEON intrinsics. Round constants are 0 as the function mimics Intel AVX2 implementation,
// which applies in order: ShiftRows, SubBytes, MixColumns, AddRoundKey.
// vaeseq_u8 applies SubBytes, ShiftRows, AddRoundKey.
static void arm_aes128_encrypt(const uint8x16_t rkeys[11], uint8x16_t nv, unsigned char *out) {

    uint8x16_t temp = vaeseq_u8(nv, rkeys[0]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[1]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[2]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[3]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[4]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[5]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[6]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[7]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[8]);

    temp = vaesmcq_u8(temp);
    temp = vaeseq_u8(temp, rkeys[9]);

    temp = veorq_u8(temp, rkeys[10]);
    vst1q_u8(out, temp);
}

void arm_aes128_free_schedule(void *schedule) {
    if (schedule != NULL) {
        mayo_secure_free(schedule, 11 * sizeof(uint16x8_t));
    }
}


static void arm_aes128_ctr_enc_sch(const void *schedule, uint8_t *out,
                                        size_t out_len) {
    uint8x16_t mask = {0, 1, 2, 3, 4, 5, 6, 7, 15, 14, 13, 12, 11, 10, 9, 8};
    uint8x16_t block = vdupq_n_u8(0); // Initialize block to zero
    while (out_len >= 16) {
        arm_aes128_encrypt(schedule, block, out);
        out += 16;
        out_len -= 16;
        block = vqtbl1q_u8(vreinterpretq_u8_u64(vaddq_u64(vreinterpretq_u64_u8(vqtbl1q_u8(block, mask)), (uint64x2_t) {0,1}) ), mask);
    }
    if (out_len > 0) {
        uint8_t tmp[16];
        arm_aes128_encrypt(schedule, block, tmp);
        memcpy(out, tmp, out_len);
    }
}

int AES_128_CTR_NEON(unsigned char *output, size_t outputByteLen,
                   const unsigned char *input) {
    void *schedule = NULL;
    arm_aes128_load_schedule(input, &schedule);
    arm_aes128_ctr_enc_sch(schedule, output, outputByteLen);
    arm_aes128_free_schedule(schedule);
    return (int)outputByteLen;
}

#endif
