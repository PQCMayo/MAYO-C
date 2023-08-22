// SPDX-License-Identifier: Apache-2.0

/**
 * Test case for sample_solution failing in crypto_sign.
 */
#include <api.h>
#include <mem.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <rng.h>

static int test_sample_solution(void) {

#if CRYPTO_BYTES == 321
    unsigned char entropy_input[48] = { 0x69,0x9e,0xce,0x02,0x03,0x00,0x00,0x80,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 };
#else
    #error "variant not supported"
#endif

    unsigned long long msglen = 32;
    unsigned long long smlen = CRYPTO_BYTES + msglen;

    unsigned char *pk  = calloc(CRYPTO_PUBLICKEYBYTES, 1);
    unsigned char *sk  = calloc(CRYPTO_SECRETKEYBYTES, 1);

    unsigned char *sig = calloc(smlen, 1);

    unsigned char msg[32] = { 0xe };
    unsigned char msgOpen[32] = { 0 };

    int res = 0;
    randombytes_init(entropy_input, NULL, 256);

    res = crypto_sign_keypair(pk, sk);
    if (res) {
        printf("crypto_sign_keypair failed\n");
        goto err;
    }

    res = crypto_sign(sig, &smlen, msg, msglen, sk);

    if (res) {
        printf("crypto_sign failed\n");
        goto err;
    }

    res = crypto_sign_open(msgOpen, &msglen, sig, smlen, pk);
    if (res || memcmp(msg, msgOpen, msglen)) {
        printf("crypto_sign_open failed\n");
        res = -1;
        goto err;
    }

err:
    free(pk);
    mayo_secure_free(sk, CRYPTO_SECRETKEYBYTES);
    free(sig);
    return res;
}

int main(int argc, char *argv[]) {
    return test_sample_solution();
}
