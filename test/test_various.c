// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <randombytes.h>
#include <time.h>
#include <mayo.h>
#include <stdalign.h>
#include <arithmetic.h>
#include <fips202.h>
#include <aes_ctr.h>
#include <simple_arithmetic.h>

#include <mayo.c>

#ifdef ENABLE_CT_TESTING
#include <valgrind/memcheck.h>
#endif
#include <debug_bench_tools.h>

int test_mayo(const mayo_params_t *p) {
    unsigned char _pk[CPK_BYTES_MAX + 1] = {0};  
    unsigned char _sk[CSK_BYTES_MAX + 1] = {0};
    unsigned char _sig[SIG_BYTES_MAX + 32 + 1] = {0};
    unsigned char _msg[32+1] = { 0 };

    // Enforce unaligned memory addresses
    unsigned char *pk  = (unsigned char *) ((uintptr_t)_pk | (uintptr_t)1);
    unsigned char *sk  = (unsigned char *) ((uintptr_t)_sk | (uintptr_t)1);
    unsigned char *sig = (unsigned char *) ((uintptr_t)_sig | (uintptr_t)1);
    unsigned char *msg = (unsigned char *) ((uintptr_t)_msg | (uintptr_t)1);

    for (int i = 0; i < 32; i++)
    {
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
    VALGRIND_MAKE_MEM_DEFINED(pk, p->cpk_bytes);
#endif

    size_t smlen = PARAM_sig_bytes(p) + 32;

    res = mayo_sign(p, sig, &smlen, msg, 32, sk);
    if (res != MAYO_OK) {
        res = -1;
        printf("sign failed!\n");
        goto err;
    }

    /*
    printf("pk: ");
    print_hex(pk, p->cpk_bytes, 40);
    printf("sk: ");
    print_hex(sk, p->csk_bytes, 40);
    printf("sm: ");
    print_hex(sig, smlen, 40); */

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

static void print_hex(const unsigned char *hex, int len, int row_len) {
    int ctr = 0;
    printf("%4d:", ctr);
    for (int i = 0; i < len;  ++i) {
        if(i % row_len == 0 && i > 0){
            ctr ++;
            printf("\n%4d:", ctr);
        }
        printf("%2x", hex[i]);
    }
    printf("\n");
}

// Evaluate public map at oil vector as sanity check
int test_eval_oil(const mayo_params_t *p) {

    int ret = MAYO_OK;
    unsigned char seed_sk[CSK_BYTES_MAX];
    unsigned char S[PK_SEED_BYTES_MAX + O_BYTES_MAX];
    uint64_t P[P1_LIMBS_MAX + P2_LIMBS_MAX];
    uint64_t P3[O_MAX * O_MAX * M_VEC_LIMBS_MAX] = {0};

    unsigned char cpk[CPK_BYTES_MAX];

    unsigned char *seed_pk;
    unsigned char O[V_MAX*O_MAX];

    const int m_vec_limbs = PARAM_m_vec_limbs(p);
    (void) m_vec_limbs;
    const int param_m = PARAM_m(p);
    const int param_k = PARAM_k(p);
    const int param_n = PARAM_n(p);
    const int param_v = PARAM_v(p);
    const int param_o = PARAM_o(p);
    const int param_O_bytes = PARAM_O_bytes(p);
    const int param_pk_seed_bytes = PARAM_pk_seed_bytes(p);
    const int param_sk_seed_bytes = PARAM_sk_seed_bytes(p);

    // seed_sk $←- B^(sk_seed bytes)
    randombytes(seed_sk, param_sk_seed_bytes);

    // S ← shake256(seedsk, pk seed bytes + O bytes)
    shake256(S, param_pk_seed_bytes + param_O_bytes, seed_sk,
             param_sk_seed_bytes);
    // seed_pk ← s[0 : pk_seed_bytes]
    seed_pk = S;

    // o ← Decode_o(s[pk_seed_bytes : pk_seed_bytes + o_bytes])
    decode(S + param_pk_seed_bytes, O, param_v * param_o);

    //memset(O, 0, param_v*param_o);
    //O[0] = 1;

    expand_P1_P2(p, P, seed_pk);

    uint64_t *P1 = P;
    uint64_t *P2 = P + PARAM_P1_limbs(p);

    //memset(P1 + 8*5, 0, sizeof(uint64_t[PARAM_P1_limbs(p)-5]));
    //memset(P2, 0, sizeof(uint64_t[PARAM_P2_limbs(p)]));

    uint64_t P2_copy[P2_LIMBS_MAX];
    memcpy(P2_copy, P2, PARAM_P2_limbs(p)*sizeof(uint64_t));

    compute_P3(p, P1, P2_copy, O, P3);

    // store seed_pk in cpk
    memcpy(cpk, seed_pk, param_pk_seed_bytes);

    uint64_t P3_upper[P3_LIMBS_MAX];
    memset(P3_upper, 1, PARAM_P3_limbs(p)*sizeof(int64_t));

    // compute Upper(P3) and store in cpk
    m_upper(p, P3, P3_upper, param_o);

    // sample s from O^k
    unsigned char s[N_MAX*K_MAX] = {0};
    for (int a = 0; a < param_k; a++)
    {    
        for (int j = 0; j < param_o; j++)
        {
            s[param_n*a + param_v + j] = rand() % 16;
            for(int i = 0; i < param_v; i++)
            {
                s[param_n*a + i] ^= mul_f(O[i*param_o + j], s[param_n*a + param_v + j]);
            }
        }
    }

    // evaluate Public map at s
    unsigned char y[2*M_MAX] = {0};
    eval_public_map(p, s, P1,P2,P3_upper, y);

    for (int i = 0; i < param_m; i++)
    {
        if(y[i] != 0){
            ret = 1;
        }
    }

    if(ret){
        printf("sanity check failed: evaluation at a random oil vector should be zero:\n");

        printf("s:\n");
        print_hex(s, param_k*param_n, param_n);

        printf("P(s):\n");
        print_hex(y, param_m, param_m);
    }

    printf("oil evaluation check result: %d (0 = OK)\n", ret);

    return ret;
}


// Check if public map satisfies P(a*x) = a^2 P(x)
int test_eval_quad(const mayo_params_t *p) {

    int ret = MAYO_OK;
    unsigned char seed_sk[CSK_BYTES_MAX];
    unsigned char S[PK_SEED_BYTES_MAX + O_BYTES_MAX];
    uint64_t P[P1_LIMBS_MAX + P2_LIMBS_MAX];
    uint64_t P3[O_MAX * O_MAX * M_VEC_LIMBS_MAX] = {0};

    unsigned char *seed_pk;
    unsigned char O[V_MAX*O_MAX];

    const int m_vec_limbs = PARAM_m_vec_limbs(p);
    (void) m_vec_limbs;

    const int param_m = PARAM_m(p);
    const int param_k = PARAM_k(p);
    const int param_n = PARAM_n(p);
    const int param_v = PARAM_v(p);
    const int param_o = PARAM_o(p);
    const int param_O_bytes = PARAM_O_bytes(p);
    const int param_pk_seed_bytes = PARAM_pk_seed_bytes(p);
    const int param_sk_seed_bytes = PARAM_sk_seed_bytes(p);

    // seed_sk $←- B^(sk_seed bytes)
    randombytes(seed_sk, param_sk_seed_bytes);

    // S ← shake256(seedsk, pk seed bytes + O bytes)
    shake256(S, param_pk_seed_bytes + param_O_bytes, seed_sk,
             param_sk_seed_bytes);
    // seed_pk ← s[0 : pk_seed_bytes]
    seed_pk = S;

    // o ← Decode_o(s[pk_seed_bytes : pk_seed_bytes + o_bytes])
    decode(S + param_pk_seed_bytes, O, param_v * param_o);

    //memset(O, 0, param_v*param_o);
    //O[1] = 1;

    expand_P1_P2(p, P, seed_pk);

    uint64_t *P1 = P;
    uint64_t *P2 = P + PARAM_P1_limbs(p);

    uint64_t P2_copy[P2_LIMBS_MAX];
    memcpy(P2_copy, P2, PARAM_P2_limbs(p)*sizeof(uint64_t));
    compute_P3(p, P1, P2_copy, O, P3);

    uint64_t P3_upper[P3_LIMBS_MAX];
    memset(P3_upper, 1, PARAM_P3_limbs(p)*sizeof(int64_t));

    // compute Upper(P3) and store in cpk
    m_upper(p, P3, P3_upper, param_o);

    // sample s 
    unsigned char s[N_MAX*K_MAX] = {0};
    for (int a = 0; a < param_k*param_n; a++)
    {    
        s[a] = rand() % 16;
    }

    unsigned char y[2*M_MAX] = {0};
    eval_public_map(p, s, P1,P2,P3_upper, y);

    for (int i = 0; i < 16; i++)
    {
        unsigned char ss[N_MAX*K_MAX] = {0};
        unsigned char yy[2*M_MAX] = {0};
        unsigned char yyy[2*M_MAX] = {0};
        
        for (int j = 0; j < param_k*param_n; j++)
        {
            ss[j] = mul_f(s[j], i);
        }

        eval_public_map(p, ss, P1,P2,P3_upper, yy);
        unsigned char i_sq = mul_f(i,i);

        for (int j = 0; j < param_m; j++)
        {
            ss[j] = mul_f(s[j], i);
            yyy[j] = mul_f(i_sq, y[j]);
            if ( yyy[j] != yy[j]){
                ret = MAYO_ERR;
            }
        }

        if (ret != 0){
            printf("%2d^2*P(s):", i);
            print_hex(yy, param_m, param_m);
            printf("P(%2d*s):  ", i);
            print_hex(yyy, param_m, param_m);

            for (int j = 0; j < param_m; j++)
            {
                yyy[j] ^= yy[j];
            }

            printf("difference:");
            print_hex(yyy, param_m, param_m);
            
        }
    }

    printf("quad evaluation check result: %d (0 = OK)\n", ret);

    return ret;
}

// Compute A in two ways and compare result
int test_A(const mayo_params_t *p) {

    int ret = MAYO_OK;
    unsigned char seed_sk[CSK_BYTES_MAX];
    unsigned char S[PK_SEED_BYTES_MAX + O_BYTES_MAX];
    uint64_t P[P1_LIMBS_MAX + P2_LIMBS_MAX];
    uint64_t P3[O_MAX * O_MAX * M_VEC_LIMBS_MAX] = {0};

    unsigned char *seed_pk;
    unsigned char O[V_MAX*O_MAX];

    const int m_vec_limbs = PARAM_m_vec_limbs(p);
    (void) m_vec_limbs;

    const int param_m = PARAM_m(p);
    const int param_k = PARAM_k(p);
    const int param_n = PARAM_n(p);
    const int param_v = PARAM_v(p);
    const int param_o = PARAM_o(p);
    const int param_O_bytes = PARAM_O_bytes(p);
    const int param_pk_seed_bytes = PARAM_pk_seed_bytes(p);
    const int param_sk_seed_bytes = PARAM_sk_seed_bytes(p);
    uint64_t Mtmp[K_MAX * O_MAX * M_VEC_LIMBS_MAX] = {0};

    // seed_sk $←- B^(sk_seed bytes)
    randombytes(seed_sk, param_sk_seed_bytes);

    // S ← shake256(seedsk, pk seed bytes + O bytes)
    shake256(S, param_pk_seed_bytes + param_O_bytes, seed_sk,
             param_sk_seed_bytes);
    // seed_pk ← s[0 : pk_seed_bytes]
    seed_pk = S;

    // o ← Decode_o(s[pk_seed_bytes : pk_seed_bytes + o_bytes])
    decode(S + param_pk_seed_bytes, O, param_v * param_o);

    expand_P1_P2(p, P, seed_pk);

    uint64_t *P1 = P;
    uint64_t *P2 = P + PARAM_P1_limbs(p);

    uint64_t P2_copy[P2_LIMBS_MAX];
    memcpy(P2_copy, P2, PARAM_P2_limbs(p)*sizeof(uint64_t));
    compute_P3(p, P1, P2_copy, O, P3);

    uint64_t P3_upper[P3_LIMBS_MAX];
    memset(P3_upper, 1, PARAM_P3_limbs(p)*sizeof(int64_t));

    // compute Upper(P3) and store in cpk
    m_upper(p, P3, P3_upper, param_o);

    // compute L_i = (P1 + P1^t)*O + P2
    uint64_t L[P2_LIMBS_MAX] = {0};
    memcpy(L, P2, PARAM_P2_limbs(p)*sizeof(uint64_t));
    P1P1t_times_O(p, P1, O, L);

    // sample vinegar
    unsigned char v[V_MAX * K_MAX] = {0};
    for (int i = 0; i < param_v*param_k; i++)
    {
        v[i] = rand() % 16;
    }

    unsigned char A1[((M_MAX+7)/8*8) * (K_MAX * O_MAX + 1)] = {0};
    unsigned char A2[M_MAX * (K_MAX * O_MAX + 1)] = {0};

    uint64_t VPV[ PARAM_m_vec_limbs(p)*V_MAX*V_MAX];
    compute_M_and_VPV(p, v, L, P1, Mtmp, VPV);

    compute_A(p, Mtmp, A1);

    unsigned char s_base[N_MAX*K_MAX] = {0};
    for (int i = 0; i < param_k; i++)
    {
        for (int j = 0; j < param_v; j++)
        {
            s_base[i*param_n + j] = v[i*param_v + j];
        }        
    }
    
    unsigned char y_base[M_MAX] = {0};
    eval_public_map(p, s_base, P1, P2, P3_upper, y_base);

    for (int i = 0; i < param_k; i++)
    {
        for (int j = 0; j < param_o; j++)
        {
            unsigned char s[N_MAX*K_MAX] = {0};
            memcpy(s, s_base, param_n*param_k);

            for (int k = 0; k < param_v; k++)
            {
                s[k + i * param_n] ^= O[k*param_o + j];
            }
            s[param_v + i*param_n + j] ^= 1;

            unsigned char y[M_MAX] = {0};
            eval_public_map(p, s, P1, P2, P3_upper, y);
            
            for (int k = 0; k < param_m; k++)
            {
                y[k] ^= y_base[k];
            }

            for (int k = 0; k < param_m; k++)
            {
                A2[k*(param_k*param_o+1) + i*param_o + j] = y[k];
            }         
        }
    }
    
    for (int i = 0; i < param_m*(param_k*param_o+1); i++)
    {
        if (A1[i] != A2[i]){
            ret = 1;

            printf("A differs in row %d, col %d\n", i / (param_k*param_o + 1), i % (param_k*param_o + 1) );
            break;
        }
    }

    if(ret == 1){

    printf("A computed normally:\n");
        for (int i = 0; i < param_m; i++)
        {
            printf("%3d: ", i);
            for (int j = 0; j < param_k*param_o + 1; j++)
            {
                if (j % param_o == 0 && j > 0){
                    printf("|");
                }
                printf("%1x", A1[i*(param_k*param_o+1)+j]);
            }
            printf("\n");
        }

        printf("A computed via evaluations:\n");
        for (int i = 0; i < param_m; i++)
        {
            printf("%3d: ", i);
            for (int j = 0; j < param_k*param_o + 1; j++)
            {
                if (j % param_o == 0 && j > 0){
                    printf("|");
                }
                printf("%1x", A2[i*(param_k*param_o+1)+j]);
            }
            printf("\n");
        }
    }

    unsigned char s[N_MAX*K_MAX] = {0};
    memcpy(s, s_base, param_n*param_k);

    // sample random s = v + o
    for (int i = 0; i < param_k; i++)
    {
        for (int j = 0; j < param_o; j++)
        {
            s[i*param_n + param_v + j] = rand() % 16;

            for (int k = 0; k < param_v; k++)
            {
                s[k + i * param_n] ^= mul_f(s[i*param_n + param_v + j], O[k*param_o + j]);
            }
        }
    }

    unsigned char y[M_MAX] = {0};
    unsigned char y1[M_MAX] = {0};
    unsigned char y2[M_MAX] = {0};
    eval_public_map(p, s, P1, P2, P3_upper, y);

    // compute y_base + A*o
    memcpy(y1, y_base, param_m);
    memcpy(y2, y_base, param_m);

    for (int i = 0; i < param_k; i++)
    {
        for (int j = 0; j < param_o; j++)
        {
            for (int k = 0; k < param_m; k++)
            {
                y1[k] ^= mul_f( s[i*param_n + param_v + j], A1[k*(param_k*param_o+1)+i*param_o+j] );
                y2[k] ^= mul_f( s[i*param_n + param_v + j], A2[k*(param_k*param_o+1)+i*param_o+j] );
            }
        }
    }    

    if (memcmp(y,y1, param_m) != 0){
        printf("A1 is wrong!\n");
        ret = MAYO_ERR;
    }

    if (memcmp(y,y2, param_m) != 0){
        printf("A2 is wrong!\n");
        ret = MAYO_ERR;
    }

    unsigned char t[M_MAX] = {0};
    unsigned char tt[M_MAX] = {0};
    unsigned char x[V_MAX*K_MAX] = {0};
    unsigned char r[V_MAX*K_MAX] = {0};

    for (int i = 0; i < param_k*param_o; i++)
    {
        r[i] = rand() % 16;
    }

    for (int i = 0; i < param_m; i++)
    {
        t[i]  = rand() % 16;
        tt[i] = t[i];
    }

    for (int i = 0; i < param_m; i++)
    {
        tt[i] ^= y_base[i];
    }
    
    sample_solution(p, A1, tt, r, x, param_k, param_o, param_m, param_k*param_o + 1);

    memcpy(s, s_base, param_n*param_k);

    unsigned char *vi;
    unsigned char Ox[V_MAX]  = {0};       

    for (int i = 0; i <= param_k - 1; ++i) {
        vi = s_base + i * param_n;
        mat_mul(O, x + i * param_o, Ox, param_o, param_v, 1);
        mat_add(vi, Ox, s + i * param_n, param_v, 1);
        memcpy(s + i * param_n + param_v, x + i * param_o, param_o);
    }

    unsigned char eval[M_MAX] = {0};
    eval_public_map(p, s, P1, P2, P3_upper, eval);


    if (memcmp(t, eval, param_m) != 0){
        ret = MAYO_ERR;
        
        printf("target:\n");
        print_hex(t, param_m, param_m);
        
        printf("eval: \n");
        print_hex(eval, param_m, param_m);
    }

    printf("A check result: %d (0 = OK)\n", ret);
    return ret;
}


// test sample_solution
int test_sample_sol(const mayo_params_t *p) {

    int ret = MAYO_OK;

    const int param_m = PARAM_m(p);
    const int param_k = PARAM_k(p);
    const int param_o = PARAM_o(p);

#define M_MAX_ROUND_UP (((M_MAX +7)/8)*8)

    unsigned char A[M_MAX_ROUND_UP * (K_MAX * O_MAX + 1)] = {0};
    unsigned char A1[M_MAX_ROUND_UP * (K_MAX * O_MAX + 1)] = {0};
    unsigned char y[M_MAX] = {0};
    unsigned char x[K_MAX * N_MAX] = {0};
    unsigned char r[K_MAX * O_MAX + 1] = {0};

    const int cols = param_k*param_o + 1;

    for (int i = 0; i < param_m*(param_k*param_o + 1); i++)
    {
        A[i] = rand() % 16;
        A1[i] = A[i];
    }
    
    for (int i = 0; i < param_m; i++)
    {
        y[i] = rand() % 16;
    }

    for (int i = 0; i < param_o*param_k; i++)
    {
        r[i] = rand() % 16;
    }

    //sample random A 

    sample_solution(p, A, y, r, x, param_k, param_o, param_m, cols);

    unsigned char yy[M_MAX] = {0};

    // compute A*x

    for (int i = 0; i < param_k*param_o; i++)
    {
        for (int k = 0; k < param_m; k++)
        {
            yy[k] ^= mul_f( x[i], A1[k*(param_k*param_o+1)+i] );
        }
    }    

    if (memcmp(y,yy, param_m) != 0){
        printf("sample solution is wrong!\n");
        ret = MAYO_ERR;

        printf("y and A*x:\n");
        print_hex(y, param_m, param_m);
        print_hex(yy, param_m, param_m);
    }

    printf("Sample solution check result: %d (0 = OK)\n", ret);
    return ret;
}

int main(int argc, char *argv[]) {
    int rc = 0;

    srand(time(NULL));

    const mayo_params_t *p;

#ifdef ENABLE_PARAMS_DYNAMIC
    if (!strcmp(argv[1], "MAYO-1")) {
        p = &MAYO_1;
    } else if (!strcmp(argv[1], "MAYO-2")) {
        p = &MAYO_2;
    } else if (!strcmp(argv[1], "MAYO-3")) {
        p = &MAYO_3;
    } else if (!strcmp(argv[1], "MAYO-5")) {
        p = &MAYO_5;
    } else {
        printf("unknown parameter set\n");
        return MAYO_ERR;
    }
#else
    p = NULL;
#endif

    rc ^= test_mayo(p);
    rc ^= test_eval_oil(p);
    rc ^= test_eval_quad(p);
    rc ^= test_A(p);
    rc ^= test_sample_sol(p);

    if (rc != MAYO_OK) {
        printf("test failed for %s\n", argv[1]);
    }
    return rc;
}