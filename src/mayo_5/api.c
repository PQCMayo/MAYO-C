// SPDX-License-Identifier: Apache-2.0

#include <api.h>
#include <mayo.h>

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
    return mayo_keypair(&MAYO_5, pk, sk);
}

#ifndef PQM4
int
crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk) {
    return mayo_sign(&MAYO_5, sm, smlen, m, mlen, sk);
}

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk) {
    return mayo_open(&MAYO_5, m, mlen, sm, smlen, pk);
}


#else
int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk) {

    unsigned long long smlen_ll;
    int rc = mayo_sign(&MAYO_5, sm, &smlen_ll, m, mlen, sk);
    *smlen = smlen_ll;
    return rc;
}

int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk) {
    unsigned long long mlen_ll;
    int rc = mayo_open(&MAYO_5, m, &mlen_ll, sm, smlen, pk);
    *mlen = mlen_ll;
    return rc;
}
#endif