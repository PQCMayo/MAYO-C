// SPDX-License-Identifier: Apache-2.0

/**
 * An example to demonstrate how to use the MAYO API.
 */

#include <mayo.h>
#include <mem.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/**
 * Example for MAYO variant:
 * - mayo_keypair
 * - mayo_expand_pk
 * - mayo_expand_sk
 * - mayo_sign
 * - mayo_verify
 * - mayo_open
 * 
 * @param[in] p MAYO parameter
 * @return int return code
 */
static int example_mayo(const mayo_params_t* p) {

    unsigned long long msglen = 32;
    unsigned long long smlen = p->sig_bytes + msglen;

    unsigned char *pk  = calloc(p->cpk_bytes, 1);
    unsigned char *sk  = calloc(p->csk_bytes, 1);

    unsigned char *epk = calloc(p->epk_bytes, 1);
    sk_t *esk = calloc(sizeof(sk_t), 1);

    unsigned char *sig = calloc(p->sig_bytes + msglen, 1);

    unsigned char msg[32] = { 0xe };
    unsigned char msg2[32] = { 0 };

    printf("Example with %s\n", p->name);

    printf("mayo_keypair -> ");
    int res = mayo_keypair(p, pk, sk);
    if (res != MAYO_OK) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("mayo_expand_sk -> ");
    res = mayo_expand_sk(p, sk, esk);
    if (res != MAYO_OK) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("mayo_expand_pk -> ");
    res = mayo_expand_pk(p, pk, epk);
    if (res != MAYO_OK) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("mayo_sign -> ");
    res = mayo_sign(p, sig, &smlen, msg, msglen, sk);
    if (res != MAYO_OK) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("mayo_open (with correct signature) -> ");
    res = mayo_open(p, msg2, &msglen, sig, smlen, pk);
    if (res != MAYO_OK || memcmp(msg, msg2, msglen)) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = MAYO_OK;
        printf("OK\n");
    }

    printf("mayo_verify (with correct signature) -> ");
    res = mayo_verify(p, msg, msglen, sig, pk);
    if (res != MAYO_OK) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = MAYO_OK;
        printf("OK\n");
    }

    printf("mayo_open (with altered signature) -> ");
    sig[0] = ~sig[0];
    memset(msg2, 0, msglen);
    res = mayo_open(p, msg2, &msglen, sig, smlen, pk);
    if (res != MAYO_ERR || !memcmp(msg, msg2, msglen)) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = MAYO_OK;
        printf("OK\n");
    }

    printf("mayo_verify (with altered signature) -> ");
    res = mayo_verify(p, msg, msglen, sig, pk);
    if (res == MAYO_OK) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = MAYO_OK;
        printf("OK\n");
    }

err:
    free(pk);
    free(epk);
    mayo_secure_free(sk, p->csk_bytes);
    free(esk);
    free(sig);
    return res;
}

int main(void) {
#ifdef ENABLE_PARAMS_DYNAMIC
    const mayo_params_t* params[] = {&MAYO_1, &MAYO_2, &MAYO_3, &MAYO_5};
    int paramslen = sizeof(params) / sizeof(params[0]);
    for (int i = 0; i < paramslen; ++i) {
        if (example_mayo(params[i]) != MAYO_OK) {
            return -1;
        }
    }
#else
    return example_mayo(&MAYO_VARIANT);
#endif
}
