// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <rng.h>
#include <mayo.h>
#include <stdalign.h>

#ifdef ENABLE_CT_TESTING
#include <valgrind/memcheck.h>
#endif

#ifdef ENABLE_CT_TESTING
static void print_hex(const unsigned char *hex, int len) {
    unsigned char *copy  = calloc(len, 1);
    memcpy(copy, hex, len); // make a copy that we can tell valgrind is okay to leak
    VALGRIND_MAKE_MEM_DEFINED(copy, len);

    for (int i = 0; i < len;  ++i) {
        printf("%02x", copy[i]);
    }
    printf("\n");
    free(copy);
}
#else
static void print_hex(const unsigned char *hex, int len) {
    for (int i = 0; i < len;  ++i) {
        printf("%02x", hex[i]);
    }
    printf("\n");
}
#endif


static int test_mayo(const mayo_params_t *p) {
    alignas(32) unsigned char pk[CPK_BYTES_MAX] = {0};  
    alignas(32) unsigned char sk[CSK_BYTES_MAX] = {0};
    alignas(32) unsigned char sig[5000 + 32] = {0};

    unsigned char seed[48] = { 0 };
    unsigned char msg[32] = { 0 };
    unsigned long long msglen = 32;

    randombytes_init(seed, NULL, 256);

    printf("Testing Keygen, Sign, Open: %s\n", p->name);

    int res = mayo_keypair(p, pk, sk);
    if (res != MAYO_OK) {
        res = -1;
        goto err;
    }

#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(pk, p->cpk_bytes);
#endif

    unsigned long long smlen = p->sig_bytes + 32;

    res = mayo_sign(p, sig, &smlen, msg, 32, sk);
    if (res != MAYO_OK) {
        res = -1;
        goto err;
    }

    printf("pk: ");
    print_hex(pk, p->cpk_bytes);
    printf("sk: ");
    print_hex(sk, p->csk_bytes);
    printf("sm: ");
    print_hex(sig, smlen);

#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(sig, smlen);
#endif

    res = mayo_open(p, msg, &msglen, sig, smlen, pk);
    if (res != MAYO_OK) {
        res = -1;
        goto err;
    }

    sig[0] = ~sig[0];
    res = mayo_open(p, msg, &msglen, sig, smlen, pk);
    if (res != MAYO_ERR) {
        res = -1;
        goto err;
    } else {
        res = MAYO_OK;
    }

err:
    return res;
}

int main(int argc, char *argv[]) {
    int rc = 0;

#ifdef ENABLE_PARAMS_DYNAMIC
    if (!strcmp(argv[1], "MAYO_1")) {
        rc = test_mayo(&MAYO_1);
    } else if (!strcmp(argv[1], "MAYO_2")) {
        rc = test_mayo(&MAYO_2);
    } else if (!strcmp(argv[1], "MAYO_3")) {
        rc = test_mayo(&MAYO_3);
    } else if (!strcmp(argv[1], "MAYO_5")) {
        rc = test_mayo(&MAYO_5);
    }
#else
    rc = test_mayo(&MAYO_VARIANT);
#endif

    if (rc != MAYO_OK) {
        printf("test failed for %s\n", argv[1]);
    }
    return rc;
}
