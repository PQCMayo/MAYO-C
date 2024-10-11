// SPDX-License-Identifier: Apache-2.0

#ifndef SHUFFLE_ARITHMETIC_64_H
#define SHUFFLE_ARITHMETIC_64_H

#include <arm_neon.h>
#include <stdint.h>
#include <mayo.h>
#include <arithmetic_common.h>

// P1*0 -> P1: v x v, O: v x o
static
inline void mayo_12_P1_times_O_neon(const uint64_t *_P1, uint8x16_t *O_multabs, uint64_t *_acc){

    const uint8x16_t *P1 = (uint8x16_t *) _P1;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    size_t cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[2*O_MAX] = {0};
        for (size_t c = r; c < V_MAX; c++)
        {
            uint8x16_t in_odd0 = P1[cols_used];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            cols_used ++;
            uint8x16_t in_odd1 = P1[cols_used];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;
            cols_used ++;

            for (size_t k = 0; k < O_MAX; k+=2)
            {
                temp[2*k]     ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd0);
                temp[2*k + 1] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even0);
                temp[2*k + 2] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd1);
                temp[2*k + 3] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even1);
            }
        }

        // convert to normal format and add to accumulator
        for (size_t k = 0; k < O_MAX; k+=2)
        {
            acc[(2*r*O_MAX) + 2*k]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[(2*r*O_MAX) + 2*k + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[(2*r*O_MAX) + 2*k + 2] ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[(2*r*O_MAX) + 2*k + 3] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
    }
}


static
inline void mayo_12_Ot_times_P1O_P2_neon(const uint64_t *_P1O_P2, uint8x16_t *O_multabs, uint64_t *_acc){
    const uint8x16_t *P1O_P2 = (uint8x16_t *) _P1O_P2;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );
    for (size_t c = 0; c < O_MAX; c++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[2*O_MAX] = {0};
        for (size_t r = 0; r < V_MAX; r++)
        {
            uint8x16_t in_odd0 = P1O_P2[2*r*O_MAX + 2*c];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1O_P2[2*r*O_MAX + 2*c + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;

            for (size_t k = 0; k < O_MAX; k+=2)
            {
                temp[2*k]     ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_odd0);
                temp[2*k + 1] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_even0);
                temp[2*k + 2] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_odd1);
                temp[2*k + 3] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_even1);
            }
        }

        // convert to normal format and add to accumulator
        for (size_t k = 0; k < O_MAX; k+=2)
        {
            acc[2*(k*O_MAX) + 2*c]         ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*(k*O_MAX) + 2*c + 1]     ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[2*((k+1)*O_MAX) + 2*c]     ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[2*((k+1)*O_MAX) + 2*c + 1] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
    }
}

static
inline void mayo_12_P1P1t_times_O_neon(const uint64_t *_P1, const unsigned char *O, uint64_t *_acc){

    const uint8x16_t *P1 = (uint8x16_t *) _P1;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    uint8x16_t O_multabs[O_MAX/2*V_MAX];
    mayo_O_multabs_neon(O, O_multabs);

    size_t cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[2*O_MAX] = {0};
        cols_used += 1;
        size_t pos = r;
        for (size_t c = 0; c < r; c++)
        {
            uint8x16_t in_odd0 = P1[2*pos];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1[2*pos + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;
            pos += (V_MAX -c - 1);

            for (size_t k = 0; k < O_MAX; k+=2)
            {
                temp[2*k]     ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd0);
                temp[2*k + 1] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even0);
                temp[2*k + 2] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd1);
                temp[2*k + 3] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even1);
            }
        }

        for (size_t c = r+1; c < V_MAX; c++)
        {
            uint8x16_t in_odd0 = P1[2*cols_used];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1[2*cols_used + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;
            cols_used ++;

            for (size_t k = 0; k < O_MAX; k+=2)
            {
                temp[2*k]     ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd0);
                temp[2*k + 1] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even0);
                temp[2*k + 2] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd1);
                temp[2*k + 3] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even1);
            }
        }

        for (size_t k = 0; k < O_MAX; k+=2)
        {
            acc[2*(r*O_MAX) + 2*k]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*(r*O_MAX) + 2*k + 2] ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[2*(r*O_MAX) + 2*k + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[2*(r*O_MAX) + 2*k + 3] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
    }
}



static
inline void mayo_12_Vt_times_L_neon(const uint64_t *_L, const uint8x16_t *V_multabs, uint64_t *_acc){

    const uint8x16_t *L = (uint8x16_t *) _L;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );
    size_t k;

    for (size_t c = 0; c < O_MAX; c++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*2] = {0};
        for (size_t r = 0; r < V_MAX; r++)
        {
            uint8x16_t in_odd0 = L[2*r*O_MAX + 2*c];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = L[2*r*O_MAX + 2*c + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[4*k]     ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even1);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k+1 < K_MAX; k+=2)
        {
            acc[2*(k*O_MAX) + 2*c]         ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*((k+1)*O_MAX) + 2*c]     ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[2*(k*O_MAX) + 2*c + 1]     ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[2*((k+1)*O_MAX) + 2*c + 1] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
#if K_MAX % 2 == 1
        acc[2*(k*O_MAX) + 2*c]         ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
        acc[2*(k*O_MAX) + 2*c + 1]     ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
#endif
    }
}

static
inline void mayo_12_P1_times_Vt_neon(const uint64_t *_P1, uint8x16_t *V_multabs, uint64_t *_acc){
    size_t k,c;
    const uint8x16_t *P1 = (uint8x16_t *) _P1;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    size_t cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*2] = {0};

        for (c=r; c < V_MAX; c++)
        {
            uint8x16_t in_odd0 = P1[2*cols_used];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1[2*cols_used + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;
            cols_used ++;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[4*k]     ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_even1);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k + 1 < K_MAX; k+=2)
        {
            acc[2*(r*K_MAX) + 2*k]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*(r*K_MAX) + 2*k + 2] ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[2*(r*K_MAX) + 2*k + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[2*(r*K_MAX) + 2*k + 3] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
#if K_MAX % 2 == 1
        acc[2*(r*K_MAX) + 2*k]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
        acc[2*(r*K_MAX) + 2*k + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
#endif
    }
}

static
inline void mayo_12_Vt_times_Pv_neon(const uint64_t *_Pv, const uint8x16_t *V_multabs, uint64_t *_acc){

    const uint8x16_t *Pv = (uint8x16_t *) _Pv;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  =  vdupq_n_u8( 0xf );
    size_t k;

    for (size_t c = 0; c < K_MAX; c++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*2] = {0};
        for (size_t r = 0; r < V_MAX; r++)
        {
            uint8x16_t in_odd0 = Pv[2*r*K_MAX + 2*c];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = Pv[2*r*K_MAX + 2*c + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[4*k]     ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even1);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k+1 < K_MAX; k+=2)
        {
            acc[2*(k*K_MAX) + 2*c]         ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*((k+1)*K_MAX) + 2*c]     ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[2*(k*K_MAX) + 2*c + 1]     ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[2*((k+1)*K_MAX) + 2*c + 1] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
#if K_MAX % 2 == 1
        acc[2*(k*K_MAX) + 2*c]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
        acc[2*(k*K_MAX) + 2*c + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
#endif
    }
}

// P2*S2 -> P2: v x o, S2: o x k
static
inline void mayo_12_P1_times_S1_plus_P2_times_S2_neon(const uint64_t *_P1, const uint64_t *_P2, uint8x16_t *S1_multabs, uint8x16_t *S2_multabs, uint64_t *_acc){
    size_t k,c;
    const uint8x16_t *P1 = (uint8x16_t *) _P1;
    const uint8x16_t *P2 = (uint8x16_t *) _P2;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    size_t P1_cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*2] = {0};


        // P1 * S1
        for (c=r; c < V_MAX; c++)
        {
            uint8x16_t in_odd0 = P1[2*P1_cols_used];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1[2*P1_cols_used + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;
            P1_cols_used ++;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[4*k]     ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_even1);
            }
        }

        // P2 * S2
        for (c=0; c < O_MAX; c++)
        {
            uint8x16_t in_odd0 = P2[2*r*O_MAX + 2*c];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P2[2*r*O_MAX + 2*c + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[4*k]     ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even1);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k + 1 < K_MAX; k+=2)
        {
            acc[2*(r*K_MAX) + 2*k]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*(r*K_MAX) + 2*k + 2] ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[2*(r*K_MAX) + 2*k + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[2*(r*K_MAX) + 2*k + 3] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
#if K_MAX % 2 == 1
        acc[2*(r*K_MAX) + 2*k]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
        acc[2*(r*K_MAX) + 2*k + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
#endif
    }
}


// P3*S2 -> P3: o x o, S2: o x k // P3 upper triangular
static
inline void mayo_12_P3_times_S2_neon(const uint64_t *_P3, uint8x16_t *S2_multabs, uint64_t *_acc){
    size_t k,c;
    const uint8x16_t *P3 = (uint8x16_t *) _P3;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    size_t cols_used = 0;
    for (size_t r = 0; r < O_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*2] = {0};

        for (c=r; c < O_MAX; c++)
        {
            uint8x16_t in_odd0 = P3[2*cols_used];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P3[2*cols_used + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;
            cols_used ++;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[4*k]     ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even1);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k + 1 < K_MAX; k+=2)
        {
            acc[2*(r*K_MAX) + 2*k]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*(r*K_MAX) + 2*k + 2] ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[2*(r*K_MAX) + 2*k + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[2*(r*K_MAX) + 2*k + 3] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
#if K_MAX % 2 == 1
            acc[2*(r*K_MAX) + 2*k]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*(r*K_MAX) + 2*k + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
#endif
    }
}


static
inline void mayo_12_S1t_times_PS1_neon(const uint64_t *_PS1, uint8x16_t *S1_multabs, uint64_t *_acc){
    mayo_12_Vt_times_Pv_neon(_PS1, S1_multabs, _acc);
}

static
inline void mayo_12_S2t_times_PS2_neon(const uint64_t *_PS2, uint8x16_t *S2_multabs, uint64_t *_acc){
    const uint8x16_t *PS2 = (uint8x16_t *) _PS2;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );
    size_t k;

    for (size_t c = 0; c < K_MAX; c++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*2] = {0};
        for (size_t r = 0; r < O_MAX; r++)
        {
            uint8x16_t in_odd0 = PS2[2*r*K_MAX + 2*c];
            uint8x16_t in_even0 = in_odd0 >> 4;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = PS2[2*r*K_MAX + 2*c + 1];
            uint8x16_t in_even1 = in_odd1 >> 4;
            in_odd1 &= low_nibble_mask;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[4*k]     ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_even1);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k+1 < K_MAX; k+=2)
        {
            acc[2*(k*K_MAX) + 2*c]         ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
            acc[2*((k+1)*K_MAX) + 2*c]     ^= vsriq_n_u8(temp[2*k+1], temp[2*k], 4);
            acc[2*(k*K_MAX) + 2*c + 1]     ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
            acc[2*((k+1)*K_MAX) + 2*c + 1] ^= vsriq_n_u8(temp[2*k+3], temp[2*k+2], 4);
        }
#if K_MAX % 2 == 1
        acc[2*(k*K_MAX) + 2*c]     ^= vsliq_n_u8(temp[2*k], temp[2*k+1], 4);
        acc[2*(k*K_MAX) + 2*c + 1] ^= vsliq_n_u8(temp[2*k+2], temp[2*k+3], 4);
#endif
    }
}


#undef K_OVER_2
#endif

