// SPDX-License-Identifier: Apache-2.0

/**
 * An example to demonstrate how to use MAYO with the NIST API.
 */

#include <api.h>
#include <mem.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * Example for MAYO variant:
 * - crypto_sign_keypair
 * - crypto_sign
 * - crypto_sign_open
 * 
 * @return int return code
 */
static int example_mayo(void) {

    size_t msglen = 32;
    size_t smlen = CRYPTO_BYTES + msglen;
    size_t siglen = CRYPTO_BYTES;

    unsigned char *pk  = calloc(CRYPTO_PUBLICKEYBYTES, 1);
    unsigned char *sk  = calloc(CRYPTO_SECRETKEYBYTES, 1);

    unsigned char *sig = calloc(smlen, 1);

    unsigned char msg[32] = { 0xe };
    unsigned char msg2[32] = { 0 };

    printf("Example with %s\n", CRYPTO_ALGNAME);

    printf("crypto_sign_keypair -> ");
    int res = crypto_sign_keypair(pk, sk);
    if (res) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("crypto_sign -> ");
    res = crypto_sign(sig, &smlen, msg, msglen, sk);
    if (res) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("crypto_sign_open (with correct signature) -> ");
    res = crypto_sign_open(msg2, &msglen, sig, smlen, pk);
    if (res || memcmp(msg, msg2, msglen)) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = 0;
        printf("OK\n");
    }

    printf("crypto_sign_open (with altered signature) -> ");
    sig[0] = ~sig[0];
    memset(msg2, 0, msglen);
    res = crypto_sign_open(msg2, &msglen, sig, smlen, pk);
    if (!res || !memcmp(msg, msg2, msglen)) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = 0;
        printf("OK\n");
    }
    
    printf("crypto_sign_signature -> ");
    res = crypto_sign_signature(sig, &siglen, msg, msglen, sk);
    if (res) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("crypto_sign_verify (with correct signature) -> ");
    res = crypto_sign_verify(sig, siglen, msg, msglen, pk);
    if (res) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = 0;
        printf("OK\n");
    }

    printf("crypto_sign_verify (with altered signature) -> ");
    sig[0] = ~sig[0];
    res = crypto_sign_verify(sig, siglen, msg, msglen, pk);
    if (!res) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = 0;
        printf("OK\n");
    }

err:
    free(pk);
    mayo_secure_free(sk, CRYPTO_SECRETKEYBYTES);
    free(sig);
    return res;
}

int main(void) {
    return example_mayo();
}

