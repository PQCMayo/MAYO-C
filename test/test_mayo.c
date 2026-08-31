// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <randombytes.h>
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
    unsigned char _pk[CPK_BYTES_MAX + 1] = {0};  
    unsigned char _sk[CSK_BYTES_MAX + 1] = {0};
    unsigned char _sig[SIG_BYTES_MAX + 32 + 1] = {0};
    unsigned char _msg[32+1] = { 0 };

    // Enforce unaligned memory addresses
    unsigned char *pk  = (unsigned char *) ((uintptr_t)_pk | (uintptr_t)1);
    unsigned char *sk  = (unsigned char *) ((uintptr_t)_sk | (uintptr_t)1);
    unsigned char *sig = (unsigned char *) ((uintptr_t)_sig | (uintptr_t)1);
    unsigned char *msg = (unsigned char *) ((uintptr_t)_msg | (uintptr_t)1);

    for (int i = 0; i < 32; i++) {
        msg[i] = i;
    }

    unsigned char seed[48] = { 0 };
    size_t msglen = 32;

    randombytes_init(seed, NULL, 256);

    printf("Testing Keygen, Sign, Open: %s\n", PARAM_name(p));

    int res = mayo_keypair(p, pk, sk);
    if (res != MAYO_OK) {
        res = -1;
        printf("keygen failed!\n");
        goto err;
    }

#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(pk, PARAM_cpk_bytes(p));
#endif

    size_t smlen = PARAM_sig_bytes(p) + 32;

    res = mayo_sign(p, sig, &smlen, msg, 32, sk);
    if (res != MAYO_OK) {
        res = -1;
        printf("sign failed!\n");
        goto err;
    }

    printf("pk: ");
    print_hex(pk, PARAM_cpk_bytes(p));
    printf("sk: ");
    print_hex(sk, PARAM_csk_bytes(p));
    printf("sm: ");
    print_hex(sig, smlen);

#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(sig, smlen);
#endif

    res = mayo_open(p, msg, &msglen, sig, smlen, pk);
    if (res != MAYO_OK) {
        res = -1;
        printf("verify failed!\n");
        goto err;
    }

    printf("verify success!\n");

    // Expanded-key API: expand once, sign and verify many times. Covers the
    // materialized-L path and the packed P1/P2 layout, which the KAT does not reach,
    // and catches an expanded key left modified by an earlier call. The expanded and
    // compact paths must agree both ways.
    {
        static sk_t esk;
        static pk_t epk;
        unsigned char sm2[SIG_BYTES_MAX + 32] = {0};
        unsigned char mout[32] = {0};
        size_t sm2len = 0, ml2 = 32;
        res = mayo_expand_sk(p, sk, &esk);
        if (res != MAYO_OK) {
            res = -1; printf("expand_sk failed!\n"); goto err;
        }
        res = mayo_expand_pk(p, pk, epk.p);
        if (res != MAYO_OK) {
            res = -1; printf("expand_pk failed!\n"); goto err;
        }
        for (int j = 0; j < 3; j++) {
            msg[0] = (unsigned char)(j + 1);
            res = mayo_sign_esk(p, sm2, &sm2len, msg, 32, &esk);
            if (res != MAYO_OK) {
                res = -1; printf("sign_esk failed!\n"); goto err;
            }
#ifdef ENABLE_CT_TESTING
            VALGRIND_MAKE_MEM_DEFINED(sm2, sm2len); // signature is public
#endif
            res = mayo_open_epk(p, mout, &ml2, sm2, sm2len, epk.p);
            if (res != MAYO_OK) {
                res = -1; printf("esk signature %d did not open with epk!\n", j); goto err;
            }
            if (ml2 != 32 || memcmp(mout, msg, 32) != 0) {
                res = -1; printf("open_epk recovered the wrong message!\n"); goto err;
            }
            res = mayo_verify_epk(p, msg, 32, sm2, epk.p);
            if (res != MAYO_OK) {
                res = -1; printf("esk signature %d did not verify with epk!\n", j); goto err;
            }
            res = mayo_open(p, mout, &ml2, sm2, sm2len, pk);
            if (res != MAYO_OK) {
                res = -1; printf("esk signature %d did not open with cpk!\n", j); goto err;
            }
        }
        // a signature made with the compact key must open with the expanded one
        res = mayo_open_epk(p, mout, &ml2, sig, smlen, epk.p);
        if (res != MAYO_OK) {
            res = -1; printf("mayo_sign signature did not open with epk!\n"); goto err;
        }

        // and a corrupted signature must be rejected by both expanded entry points
        // (msg still matches sm2 here, so only the corruption can cause a reject)
        sm2[0] = ~sm2[0];
        if (mayo_verify_epk(p, msg, 32, sm2, epk.p) != MAYO_ERR ||
            mayo_open_epk(p, mout, &ml2, sm2, sm2len, epk.p) != MAYO_ERR) {
            res = -1; printf("epk accepted a corrupted signature!\n"); goto err;
        }
        msg[0] = 0;
        res = MAYO_OK;
        printf("expand_sk/expand_pk + sign_esk + open_epk: verified\n");
    }

    sig[0] = ~sig[0];
    res = mayo_open(p, msg, &msglen, sig, smlen, pk);
    if (res != MAYO_ERR) {
        res = -1;
        printf("wrong signature still verified!\n");
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
    if (!strcmp(argv[1], "MAYO-1")) {
        rc = test_mayo(&MAYO_1);
    } else if (!strcmp(argv[1], "MAYO-2")) {
        rc = test_mayo(&MAYO_2);
    } else if (!strcmp(argv[1], "MAYO-3")) {
        rc = test_mayo(&MAYO_3);
    } else if (!strcmp(argv[1], "MAYO-5")) {
        rc = test_mayo(&MAYO_5);
    } else {
        printf("unknown parameter set\n");
        return MAYO_ERR;
    }
#else
    rc = test_mayo(NULL);
#endif

    if (rc != MAYO_OK) {
        printf("test failed for %s\n", argv[1]);
    }
    return rc;
}

