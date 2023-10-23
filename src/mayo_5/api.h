// SPDX-License-Identifier: Apache-2.0

#ifndef api_h
#define api_h

#define CRYPTO_SECRETKEYBYTES 40
#define CRYPTO_PUBLICKEYBYTES 5008
#define CRYPTO_BYTES 838

#define CRYPTO_ALGNAME "MAYO_5"

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk);

#ifndef PQM4
int
crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk);

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk);

#else
#include <stddef.h>

int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk);

int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk);
#endif
#endif /* api_h */
