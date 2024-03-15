// SPDX-License-Identifier: Apache-2.0

#include <mayo.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdalign.h>
#include <rng.h>


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
      printf("  %-35s -> median: %2" PRId64 ", average: %2" PRId64 " ", name, \
      cycles_list[(runs < LIST_SIZE)? runs/2 : LIST_SIZE/2], (cycles / runs)); \
      printf("%s\n", BENCH_UNITS); \
    }

#define LIST_SIZE 10000

extern void P1_times_O(const mayo_params_t* p, const uint32_t* P1, const unsigned char* O, uint32_t* acc);
extern void mul_add_mat_trans_x_bitsliced_m_mat(int m_legs, const unsigned char *mat, const uint32_t *bs_mat, uint32_t *acc, int mat_rows, int mat_cols, int bs_mat_cols);
extern void mul_add_mat_x_bitsliced_m_mat(int m_legs, const unsigned char *mat, const uint32_t *bs_mat, uint32_t *acc, int mat_rows, int mat_cols, int bs_mat_cols);
extern void P1_times_Vt(const mayo_params_t* p, const uint32_t* P1, const unsigned char* V, uint32_t* acc);
extern void bitsliced_m_calculate_PS_SPS(const uint32_t *bitsliced_P1, const uint32_t *bitsliced_P2, const uint32_t *bitsliced_P3, const unsigned char *S,
                              const int m, const int v, const int o, const int k, uint32_t *bitsliced_SPS);

static int bench_sig(const mayo_params_t *p, int runs, int csv) {

    int rc = 0;
    int i;

    int64_t cycles, cycles1, cycles2;
    int64_t cycles_list[runs];

    // fill variables with dummy data
    alignas (32) uint32_t bitsliced_P[(P1_BYTES_MAX + P2_BYTES_MAX) / 4] = {1};
    alignas (32) uint32_t bitsliced_P3[O_MAX * O_MAX * M_MAX / 8] = {2};
    unsigned char Vdec[N_MINUS_O_MAX * K_MAX] = {3};
    uint32_t bitsliced_M[K_MAX * O_MAX * M_MAX / 8] = {4};
    unsigned char O[(N_MINUS_O_MAX)*O_MAX] = {5};
    unsigned char s[K_MAX * N_MAX] = {6};

    uint32_t *bitsliced_P1 = bitsliced_P;
    uint32_t *bitsliced_P2 = bitsliced_P + (P1_BYTES_MAX / 4);
    uint32_t *bitsliced_L = bitsliced_P + (P1_BYTES_MAX/4);
    uint32_t *bitsliced_P1O_P2 = bitsliced_P2;

    if (csv) {
        printf("%s,", p->name);
    } else {
        printf("Benchmarking %s\n", p->name);
    }

    BENCH_CODE_1(runs);
    P1_times_O(p, bitsliced_P1, O, bitsliced_P1O_P2);
    mul_add_mat_trans_x_bitsliced_m_mat(M_MAX/32, O, bitsliced_P1O_P2, bitsliced_P3,
                                        V_MAX, O_MAX, O_MAX);
    BENCH_CODE_2("Tab.3/Col 1: -O^t * (P1*O + P2)", csv);

    BENCH_CODE_1(runs);
    // compute all the v_i^t * P^(1) * v_j
    alignas (32) uint32_t bitsliced_Pv[N_MINUS_O_MAX * K_MAX * M_MAX / 8] = {0};
    alignas (32) uint32_t bitsliced_vPv[K_MAX * K_MAX * M_MAX / 8] = {0};
    // compute all the v_i^T * L matrices.
    mul_add_mat_x_bitsliced_m_mat(M_MAX / 32, Vdec, bitsliced_L, bitsliced_M,
                                      K_MAX, N_MAX - O_MAX, O_MAX);

    P1_times_Vt(p, bitsliced_P1, Vdec, bitsliced_Pv);
    mul_add_mat_x_bitsliced_m_mat(M_MAX / 32, Vdec, bitsliced_Pv,
                                      bitsliced_vPv, K_MAX, N_MAX - O_MAX,
                                      K_MAX);
    BENCH_CODE_2("Tab.3/Col 2: V*P1*V^t & V*L", csv);

    BENCH_CODE_1(runs);
    alignas (32) uint32_t bitsliced_SPS[K_MAX * K_MAX * M_MAX / 8] = {0};
    bitsliced_m_calculate_PS_SPS(bitsliced_P1, bitsliced_P2, bitsliced_P3, s, M_MAX,
                             V_MAX, O_MAX, K_MAX, bitsliced_SPS);
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
#else
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
#endif
}