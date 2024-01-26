// SPDX-License-Identifier: Apache-2.0

#include <mayo.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>


#if defined(TARGET_OS_UNIX) && (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_OTHER))
#include <time.h>
#endif
#if (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_S390X) || defined(TARGET_OTHER))
#define print_unit printf("nsec\n");
#else
#define print_unit printf("cycles\n");
#endif

static int bench_sig(const mayo_params_t *p, int runs, int csv);
static inline int64_t cpucycles(void);

int main(int argc, char *argv[]) {
    int rc = 0;

#ifdef ENABLE_PARAMS_DYNAMIC
    if (argc < 3) {
        printf("Two arguments needed\n");
        rc = 1;
        goto end;
    }
    int runs = atoi(argv[2]);
    if (!strcmp(argv[1], "MAYO_1")) {
        rc = bench_sig(&MAYO_1, runs, 0);
    } else if (!strcmp(argv[1], "MAYO_2")) {
        rc = bench_sig(&MAYO_2, runs, 0);
    } else if (!strcmp(argv[1], "MAYO_3")) {
        rc = bench_sig(&MAYO_3, runs, 0);
    } else if (!strcmp(argv[1], "MAYO_5")) {
        rc = bench_sig(&MAYO_5, runs, 0);
    }
#else
    if (argc < 2) {
        printf("One argument needed\n");
        rc = 1;
        goto end;
    }
    int runs = atoi(argv[1]);
    rc = bench_sig(&MAYO_VARIANT, runs, 0);
#endif



end:
    return rc;
}

#if (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_S390X))
#define BENCH_UNITS "nsec"
#else
#define BENCH_UNITS "cycles"
#endif

int cmpfunc (const void *a, const void *b) {
    return ( *(uint64_t *)a - * (uint64_t *)b );
}

#define BENCH_CODE_1(r) \
    cycles = 0; \
    for (i = 0; i < (r); ++i) { \
        cycles1 = cpucycles();

#define BENCH_CODE_2(name, csv) \
        cycles2 = cpucycles(); \
        if(i < LIST_SIZE) \
          cycles_list[i] = (cycles2 - cycles1);\
        cycles = cycles + (cycles2 - cycles1); \
    } \
    qsort(cycles_list, (runs < LIST_SIZE)? runs : LIST_SIZE, sizeof(uint64_t), cmpfunc);\
    if (csv) \
      printf("%2" PRId64 ",", cycles_list[(runs < LIST_SIZE)? runs/2 : LIST_SIZE/2]); \
    else { \
      printf("  %-20s-> median: %2" PRId64 ", average: %2" PRId64 " ", name, \
      cycles_list[(runs < LIST_SIZE)? runs/2 : LIST_SIZE/2], (cycles / runs)); \
      printf("%s\n", BENCH_UNITS); \
    }

#define LIST_SIZE 10000

static int bench_sig(const mayo_params_t *p, int runs, int csv) {

    int rc = 0;
    int i;

    int64_t cycles, cycles1, cycles2;
    int64_t cycles_list[10000];

    const int m_len = 32;

    unsigned char *pk  = calloc(p->cpk_bytes, 1);
    unsigned char *epk  = calloc(p->epk_bytes, 1);
    unsigned char *sk  = calloc(p->csk_bytes, 1);
    sk_t *esk  = (sk_t *)calloc(sizeof(sk_t), 1);
    unsigned char *sig = calloc(p->sig_bytes + m_len, 1);
    unsigned char *m   = calloc(m_len, 1);
    unsigned long long len = p->sig_bytes;

    if (csv) {
        printf("%s,", p->name);
    } else {
        printf("Benchmarking %s\n", p->name);
    }

    BENCH_CODE_1(runs);
    mayo_keypair(p, pk, sk);
    BENCH_CODE_2("mayo_keypair", csv);

    BENCH_CODE_1(runs);
    mayo_expand_sk(p, sk, esk);
    BENCH_CODE_2("mayo_expand_sk", csv);

    BENCH_CODE_1(runs);
    mayo_expand_pk(p, pk, epk);
    BENCH_CODE_2("mayo_expand_pk", csv);

    BENCH_CODE_1(runs);
    mayo_sign(p, sig, &len, m, m_len, sk);
    BENCH_CODE_2("mayo_sign", csv);

    len = 32;
    BENCH_CODE_1(runs);
    mayo_open(p, m, &len, sig, p->sig_bytes, pk);
    BENCH_CODE_2("mayo_verify", csv);

    if (csv) {
        printf("\n");
    }

    free(pk);
    free(sk);
    free(epk);
    free(esk);
    free(sig);
    free(m);
    return rc;
}

static inline int64_t cpucycles(void) {
#if (defined(TARGET_AMD64) || defined(TARGET_X86))
    unsigned int hi, lo;

    asm volatile ("rdtsc" : "=a" (lo), "=d"(hi));
    return ((int64_t) lo) | (((int64_t) hi) << 32);
#elif (defined(TARGET_S390X))
    uint64_t tod;
    asm volatile("stckf %0\n" : "=Q" (tod) : : "cc");
    return (tod * 1000 / 4096);
#else
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
#endif
}
