// SPDX-License-Identifier: Apache-2.0

#include <mayo.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdalign.h>

#include "m1cycles.h"

#if (defined(TARGET_OS_UNIX) && (defined(TARGET_ARM) || defined(TARGET_ARM64)) || defined(TARGET_OTHER)) \
    || (!defined(TARGET_OS_MAC) && defined(TARGET_ARM64))
#include <time.h>
#endif

static int bench_sig(const mayo_params_t *p, int runs, int csv);
static inline int64_t cpucycles(void);

int main(int argc, char *argv[]) {
    int rc = 0;

#if defined(TARGET_OS_MAC) && defined(TARGET_ARM64)
    setup_rdtsc();
#endif

#ifdef ENABLE_PARAMS_DYNAMIC
    if (argc < 3) {
        printf("Two arguments needed\n");
        rc = 1;
        goto end;
    }
    int runs = atoi(argv[2]);
    if (!strcmp(argv[1], "MAYO-1")) {
        rc = bench_sig(&MAYO_1, runs, 0);
    } else if (!strcmp(argv[1], "MAYO-2")) {
        rc = bench_sig(&MAYO_2, runs, 0);
    } else if (!strcmp(argv[1], "MAYO-3")) {
        rc = bench_sig(&MAYO_3, runs, 0);
    } else if (!strcmp(argv[1], "MAYO-5")) {
        rc = bench_sig(&MAYO_5, runs, 0);
    }
#else
    if (argc < 2) {
        printf("One argument needed\n");
        rc = 1;
        goto end;
    }
    int runs = atoi(argv[1]);
    rc = bench_sig(0, runs, 0);
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
      printf("  %-35s -> median: %2" PRId64 ", average: %2" PRId64 " ", name, \
      cycles_list[(runs < LIST_SIZE)? runs/2 : LIST_SIZE/2], (cycles / runs)); \
      printf("%s\n", BENCH_UNITS); \
    }

#define LIST_SIZE 10000

#define Ot_times_P1O_P2 MAYO_NAMESPACE(Ot_times_P1O_P2)
extern void Ot_times_P1O_P2(const mayo_params_t* p, const uint64_t* P1, const unsigned char* O, uint64_t* P1O_P2, uint64_t* P3);
#define V_times_L__V_times_P1_times_Vt MAYO_NAMESPACE(V_times_L__V_times_P1_times_Vt)
extern void V_times_L__V_times_P1_times_Vt(const mayo_params_t* p, const uint64_t* L, const unsigned char* V, uint64_t* M, const uint64_t* P1, uint64_t* Y);
#define m_calculate_PS_SPS MAYO_NAMESPACE(m_calculate_PS_SPS)
extern void m_calculate_PS_SPS(const uint64_t *P1, const uint64_t *P2, const uint64_t *P3, const unsigned char *S,
                              const int m, const int v, const int o, const int k, uint64_t *SPS);

static int bench_sig(const mayo_params_t *p, int runs, int csv) {

    int rc = 0;
    int i;

    int64_t cycles, cycles1, cycles2;
    int64_t cycles_list[runs];

    // fill variables with dummy data
    alignas (32) uint64_t P[(P1_BYTES_MAX + P2_BYTES_MAX) / 8] = {0};
    unsigned char O[(N_MINUS_O_MAX)*O_MAX] = {1};
    alignas (32) uint64_t P3[O_MAX * O_MAX * M_MAX / 16] = {2};

    alignas (32) uint64_t Y[K_MAX * K_MAX * M_MAX / 16] = {3};
    alignas (32) uint64_t Mtmp[K_MAX * O_MAX * M_MAX / 16] = {4};
    unsigned char Vdec[N_MINUS_O_MAX * K_MAX] = {5};

    alignas (64) uint64_t pkv[EPK_BYTES_MAX / 8] = {6};
    unsigned char s[K_MAX * N_MAX] = {7};
    alignas (32) uint64_t SPS[K_MAX * K_MAX * M_MAX / 16] = {8};


    if (csv) {
        printf("%s,", PARAM_name(p));
    } else {
        printf("Benchmarking %s\n", PARAM_name(p));
    }

    BENCH_CODE_1(runs);
    Ot_times_P1O_P2(0, P, O, P + (P1_BYTES_MAX / 8), P3);
    BENCH_CODE_2("Tab.3/Col 1: -O^t * (P1*O + P2)", csv);

    BENCH_CODE_1(runs);
    V_times_L__V_times_P1_times_Vt(0, P + (P1_BYTES_MAX / 8), Vdec, Mtmp, P, Y);
    BENCH_CODE_2("Tab.3/Col 2: V*P1*V^t & V*L", csv);

    BENCH_CODE_1(runs);
    m_calculate_PS_SPS(pkv, pkv + (P1_BYTES_MAX / 8), pkv + (P1_BYTES_MAX / 8) + (P2_BYTES_MAX / 8), s, M_MAX, V_MAX, O_MAX, K_MAX, SPS);
    BENCH_CODE_2("Tab.3/Col 3: S*P*S^t", csv);

    if (csv) {
        printf("\n");
    }

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
#elif (defined(TARGET_OS_MAC) && defined(TARGET_ARM64))
    return rdtsc();
#else
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
#endif
}
