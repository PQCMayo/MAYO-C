// SPDX-License-Identifier: Apache-2.0

#ifndef SHUFFLE_ARITHMETIC_64_H
#define SHUFFLE_ARITHMETIC_64_H

#include <arm_neon.h>
#include <stdint.h>
#include <mayo.h>
#include <arithmetic_common.h>

// P1*0 -> P1: v x v, O: v x o
static
inline void mayo_5_P1_times_O_neon(const uint64_t *_P1, uint8x16_t *O_multabs, uint64_t *_acc){

    const uint8x16_t *P1 = (uint8x16_t *) _P1;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    size_t cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[4*O_MAX] = {0};
        for (size_t c = r; c < V_MAX; c++)
        {
            uint8x16_t in_odd0 = P1[cols_used];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            cols_used ++;
            uint8x16_t in_odd1 = P1[cols_used];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            cols_used ++;
            uint8x16_t in_odd2 = P1[cols_used];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            cols_used ++;
            uint8x16_t in_odd3 = P1[cols_used];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;
            cols_used ++;

            for (size_t k = 0; k < O_MAX; k+=2)
            {
                temp[4*k]     ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even1);
                temp[4*k + 4] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd2);
                temp[4*k + 5] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even2);
                temp[4*k + 6] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd3);
                temp[4*k + 7] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even3);
            }
        }

        // convert to normal format and add to accumulator
        for (size_t k = 0; k < O_MAX; k+=2)
        {
            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k + 2] >> 4)) & low_nibble_mask;
            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k + 4] >> 4)) & low_nibble_mask;
            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k + 6] >> 4)) & low_nibble_mask;
            acc[(4*r*O_MAX) + 4*k]     ^= temp[4*k]     ^ (t0 << 4);
            acc[(4*r*O_MAX) + 4*k + 1] ^= temp[4*k + 2] ^ (t1 << 4);
            acc[(4*r*O_MAX) + 4*k + 2] ^= temp[4*k + 4] ^ (t2 << 4);
            acc[(4*r*O_MAX) + 4*k + 3] ^= temp[4*k + 6] ^ (t3 << 4);
            acc[(4*r*O_MAX) + 4*k + 4] ^= temp[4*k + 1] ^ t0;
            acc[(4*r*O_MAX) + 4*k + 5] ^= temp[4*k + 3] ^ t1;
            acc[(4*r*O_MAX) + 4*k + 6] ^= temp[4*k + 5] ^ t2;
            acc[(4*r*O_MAX) + 4*k + 7] ^= temp[4*k + 7] ^ t3;
        }
    }
}


static
inline void mayo_5_Ot_times_P1O_P2_neon(const uint64_t *_P1O_P2, uint8x16_t *O_multabs, uint64_t *_acc){
    const uint8x16_t *P1O_P2 = (uint8x16_t *) _P1O_P2;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );
    for (size_t c = 0; c < O_MAX; c++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[4*O_MAX] = {0};
        for (size_t r = 0; r < V_MAX; r++)
        {
            uint8x16_t in_odd0 = P1O_P2[4*r*O_MAX + 4*c];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1O_P2[4*r*O_MAX + 4*c + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = P1O_P2[4*r*O_MAX + 4*c + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = P1O_P2[4*r*O_MAX + 4*c + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;

            for (size_t k = 0; k < O_MAX; k+=2)
            {
                temp[4*k]     ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_even1);
                temp[4*k + 4] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_odd2);
                temp[4*k + 5] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_even2);
                temp[4*k + 6] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_odd3);
                temp[4*k + 7] ^= vqtbl1q_u8(O_multabs[O_MAX/2*r + k/2], in_even3);
            }
        }

        // convert to normal format and add to accumulator
        for (size_t k = 0; k < O_MAX; k+=2)
        {
            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k + 2] >> 4)) & low_nibble_mask;
            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k + 4] >> 4)) & low_nibble_mask;
            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k + 6] >> 4)) & low_nibble_mask;
            acc[4*(k*O_MAX) + 4*c]      ^= temp[4*k + 0] ^ (t0 << 4);
            acc[4*(k*O_MAX) + 4*c + 1]  ^= temp[4*k + 2] ^ (t1 << 4);
            acc[4*(k*O_MAX) + 4*c + 2]  ^= temp[4*k + 4] ^ (t2 << 4);
            acc[4*(k*O_MAX) + 4*c + 3]  ^= temp[4*k + 6] ^ (t3 << 4);
            acc[4*((k+1)*O_MAX) + 4*c]     ^= temp[4*k + 1] ^ t0;
            acc[4*((k+1)*O_MAX) + 4*c + 1] ^= temp[4*k + 3] ^ t1;
            acc[4*((k+1)*O_MAX) + 4*c + 2] ^= temp[4*k + 5] ^ t2;
            acc[4*((k+1)*O_MAX) + 4*c + 3] ^= temp[4*k + 7] ^ t3;
        }
    }
}

static
inline void mayo_5_P1P1t_times_O_neon(const uint64_t *_P1, const unsigned char *O, uint64_t *_acc){

    const uint8x16_t *P1 = (uint8x16_t *) _P1;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    uint8x16_t O_multabs[O_MAX/2*V_MAX];
    mayo_O_multabs_neon(O, O_multabs);

    size_t cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[4*O_MAX] = {0};
        cols_used += 1;
        size_t pos = r;
        for (size_t c = 0; c < r; c++)
        {
            uint8x16_t in_odd0 = P1[4*pos];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1[4*pos + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = P1[4*pos + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = P1[4*pos + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;
            pos += (V_MAX -c - 1);


            for (size_t k = 0; k < O_MAX; k+=2)
            {
                temp[4*k]     ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even1);
                temp[4*k + 4] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd2);
                temp[4*k + 5] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even2);
                temp[4*k + 6] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd3);
                temp[4*k + 7] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even3);
            }
        }

        for (size_t c = r+1; c < V_MAX; c++)
        {
            uint8x16_t in_odd0 = P1[4*cols_used];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1[4*cols_used + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = P1[4*cols_used + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = P1[4*cols_used + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;
            cols_used ++;

            for (size_t k = 0; k < O_MAX; k+=2)
            {
                temp[4*k]     ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd0);
                temp[4*k + 1] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even0);
                temp[4*k + 2] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd1);
                temp[4*k + 3] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even1);
                temp[4*k + 4] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd2);
                temp[4*k + 5] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even2);
                temp[4*k + 6] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_odd3);
                temp[4*k + 7] ^= vqtbl1q_u8(O_multabs[O_MAX/2*c + k/2], in_even3);
            }
        }

        for (size_t k = 0; k < O_MAX; k+=2)
        {
            uint8x16_t acc0 = acc[4*(r*O_MAX) + 4*k    ];
            uint8x16_t acc1 = acc[4*(r*O_MAX) + 4*k + 1];
            uint8x16_t acc2 = acc[4*(r*O_MAX) + 4*k + 2];
            uint8x16_t acc3 = acc[4*(r*O_MAX) + 4*k + 3];
            uint8x16_t acc4 = acc[4*(r*O_MAX) + 4*k + 4];
            uint8x16_t acc5 = acc[4*(r*O_MAX) + 4*k + 5];
            uint8x16_t acc6 = acc[4*(r*O_MAX) + 4*k + 6];
            uint8x16_t acc7 = acc[4*(r*O_MAX) + 4*k + 7];


            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            acc[4*(r*O_MAX) + 4*k] =     acc0 ^ temp[4*k] ^ (t0 << 4);
            acc[4*(r*O_MAX) + 4*k + 4] = acc4 ^ temp[4*k+1] ^ t0;

            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k + 2] >> 4)) & low_nibble_mask;

            acc[4*(r*O_MAX) + 4*k + 1] = acc1 ^ temp[4*k+2] ^ (t1 << 4);
            acc[4*(r*O_MAX) + 4*k + 5] = acc5 ^ temp[4*k+3] ^ t1;

            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k + 4] >> 4)) & low_nibble_mask;

            acc[4*(r*O_MAX) + 4*k + 2] = acc2 ^ temp[4*k+4] ^ (t2 << 4);
            acc[4*(r*O_MAX) + 4*k + 6] = acc6 ^ temp[4*k+5] ^ t2;

            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k + 6] >> 4)) & low_nibble_mask;

            acc[4*(r*O_MAX) + 4*k + 3] = acc3 ^ temp[4*k+6] ^ (t3 << 4);
            acc[4*(r*O_MAX) + 4*k + 7] = acc7 ^ temp[4*k+7] ^ t3;
        }
    }
}



static
inline void mayo_5_Vt_times_L_neon(const uint64_t *_L, const uint8x16_t *V_multabs, uint64_t *_acc){

    const uint8x16_t *L = (uint8x16_t *) _L;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );
    size_t k;

    for (size_t c = 0; c < O_MAX; c++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*4] = {0};
        for (size_t r = 0; r < V_MAX; r++)
        {
            uint8x16_t in_odd0 = L[4*r*O_MAX + 4*c];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = L[4*r*O_MAX + 4*c + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = L[4*r*O_MAX + 4*c + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = L[4*r*O_MAX + 4*c + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[8*k]     ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd0);
                temp[8*k + 1] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even0);
                temp[8*k + 2] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd1);
                temp[8*k + 3] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even1);
                temp[8*k + 4] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd2);
                temp[8*k + 5] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even2);
                temp[8*k + 6] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd3);
                temp[8*k + 7] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even3);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k+1 < K_MAX; k+=2)
        {
            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            acc[4*(k*O_MAX) + 4*c]     ^= temp[4*k] ^ (t0 << 4);
            acc[4*((k+1)*O_MAX) + 4*c] ^= temp[4*k+1] ^ t0;

            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k + 2] >> 4)) & low_nibble_mask;
            acc[4*(k*O_MAX) + 4*c + 1]     ^= temp[4*k+2] ^ (t1 << 4);
            acc[4*((k+1)*O_MAX) + 4*c + 1] ^= temp[4*k+3] ^ t1;

            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k + 4] >> 4)) & low_nibble_mask;
            acc[4*(k*O_MAX) + 4*c + 2]     ^= temp[4*k+4] ^ (t2 << 4);
            acc[4*((k+1)*O_MAX) + 4*c + 2] ^= temp[4*k+5] ^ t2;

            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k + 6] >> 4)) & low_nibble_mask;
            acc[4*(k*O_MAX) + 4*c + 3]     ^= temp[4*k+6] ^ (t3 << 4);
            acc[4*((k+1)*O_MAX) + 4*c + 3] ^= temp[4*k+7] ^ t3;
        }
    }
}

static
inline void mayo_5_P1_times_Vt_neon(const uint64_t *_P1, uint8x16_t *V_multabs, uint64_t *_acc){
    size_t k,c;
    const uint8x16_t *P1 = (uint8x16_t *) _P1;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    size_t cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*4] = {0};

        for (c=r; c < V_MAX; c++)
        {
            uint8x16_t in_odd0 = P1[4*cols_used];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1[4*cols_used + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = P1[4*cols_used + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = P1[4*cols_used + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;
            cols_used ++;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[8*k]     ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_odd0);
                temp[8*k + 1] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_even0);
                temp[8*k + 2] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_odd1);
                temp[8*k + 3] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_even1);
                temp[8*k + 4] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_odd2);
                temp[8*k + 5] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_even2);
                temp[8*k + 6] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_odd3);
                temp[8*k + 7] ^= vqtbl1q_u8(V_multabs[K_OVER_2*c + k], in_even3);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k + 1 < K_MAX; k+=2)
        {
            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k] ^= temp[4*k] ^ (t0 << 4);
            acc[4*(r*K_MAX) + 4*k + 4] ^= temp[4*k+1] ^ t0;

            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k+2] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 1] ^= temp[4*k+2] ^ (t1 << 4);
            acc[4*(r*K_MAX) + 4*k + 5] ^= temp[4*k+3] ^ t1;

            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k+4] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 2] ^= temp[4*k+4] ^ (t2 << 4);
            acc[4*(r*K_MAX) + 4*k + 6] ^= temp[4*k+5] ^ t2;

            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k+6] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 3] ^= temp[4*k+6] ^ (t3 << 4);
            acc[4*(r*K_MAX) + 4*k + 7] ^= temp[4*k+7] ^ t3;
        }
    }
}

static
inline void mayo_5_Vt_times_Pv_neon(const uint64_t *_Pv, const uint8x16_t *V_multabs, uint64_t *_acc){

    const uint8x16_t *Pv = (uint8x16_t *) _Pv;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  =  vdupq_n_u8( 0xf );
    size_t k;

    for (size_t c = 0; c < K_MAX; c++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*4] = {0};
        for (size_t r = 0; r < V_MAX; r++)
        {
            uint8x16_t in_odd0 = Pv[4*r*K_MAX + 4*c];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = Pv[4*r*K_MAX + 4*c + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = Pv[4*r*K_MAX + 4*c + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = Pv[4*r*K_MAX + 4*c + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[8*k]     ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd0);
                temp[8*k + 1] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even0);
                temp[8*k + 2] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd1);
                temp[8*k + 3] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even1);
                temp[8*k + 4] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd2);
                temp[8*k + 5] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even2);
                temp[8*k + 6] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_odd3);
                temp[8*k + 7] ^= vqtbl1q_u8(V_multabs[K_OVER_2*r + k], in_even3);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k+1 < K_MAX; k+=2)
        {
            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            acc[4*(k*K_MAX) + 4*c] ^= temp[4*k] ^ (t0 << 4);
            acc[4*((k+1)*K_MAX) + 4*c] ^= temp[4*k+1] ^ t0;

            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k+2] >> 4)) & low_nibble_mask;
            acc[4*(k*K_MAX) + 4*c + 1] ^= temp[4*k+2] ^ (t1 << 4);
            acc[4*((k+1)*K_MAX) + 4*c + 1] ^= temp[4*k+3] ^ t1;

            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k+4] >> 4)) & low_nibble_mask;
            acc[4*(k*K_MAX) + 4*c + 2] ^= temp[4*k+4] ^ (t2 << 4);
            acc[4*((k+1)*K_MAX) + 4*c + 2] ^= temp[4*k+5] ^ t2;

            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k+6] >> 4)) & low_nibble_mask;
            acc[4*(k*K_MAX) + 4*c + 3] ^= temp[4*k+6] ^ (t3 << 4);
            acc[4*((k+1)*K_MAX) + 4*c + 3] ^= temp[4*k+7] ^ t3;
        }
    }
}

// P2*S2 -> P2: v x o, S2: o x k
static
inline void mayo_5_P1_times_S1_plus_P2_times_S2_neon(const uint64_t *_P1, const uint64_t *_P2, uint8x16_t *S1_multabs, uint8x16_t *S2_multabs, uint64_t *_acc){
    size_t k,c;
    const uint8x16_t *P1 = (uint8x16_t *) _P1;
    const uint8x16_t *P2 = (uint8x16_t *) _P2;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    size_t P1_cols_used = 0;
    for (size_t r = 0; r < V_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*4] = {0};


        // P1 * S1
        for (c=r; c < V_MAX; c++)
        {
            uint8x16_t in_odd0 = P1[4*P1_cols_used];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P1[4*P1_cols_used + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = P1[4*P1_cols_used + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = P1[4*P1_cols_used + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;
            P1_cols_used ++;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[8*k]     ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_odd0);
                temp[8*k + 1] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_even0);
                temp[8*k + 2] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_odd1);
                temp[8*k + 3] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_even1);
                temp[8*k + 4] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_odd2);
                temp[8*k + 5] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_even2);
                temp[8*k + 6] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_odd3);
                temp[8*k + 7] ^= vqtbl1q_u8(S1_multabs[K_OVER_2*c + k], in_even3);
            }
        }

        // P2 * S2
        for (c=0; c < O_MAX; c++)
        {
            uint8x16_t in_odd0 = P2[4*r*O_MAX + 4*c];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P2[4*r*O_MAX + 4*c + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = P2[4*r*O_MAX + 4*c + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd3 = P2[4*r*O_MAX + 4*c + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[8*k]     ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd0);
                temp[8*k + 1] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even0);
                temp[8*k + 2] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd1);
                temp[8*k + 3] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even1);
                temp[8*k + 4] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd2);
                temp[8*k + 5] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even2);
                temp[8*k + 6] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd3);
                temp[8*k + 7] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even3);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k + 1 < K_MAX; k+=2)
        {
            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k] ^= temp[4*k] ^ (t0 << 4);
            acc[4*(r*K_MAX) + 4*k + 4] ^= temp[4*k+1] ^ t0;

            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k + 2] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 1] ^= temp[4*k+2] ^ (t1 << 4);
            acc[4*(r*K_MAX) + 4*k + 5] ^= temp[4*k+3] ^ t1;

            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k + 4] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 2] ^= temp[4*k+4] ^ (t2 << 4);
            acc[4*(r*K_MAX) + 4*k + 6] ^= temp[4*k+5] ^ t2;

            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k + 6] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 3] ^= temp[4*k+6] ^ (t3 << 4);
            acc[4*(r*K_MAX) + 4*k + 7] ^= temp[4*k+7] ^ t3;
        }
    }
}


// P3*S2 -> P3: o x o, S2: o x k // P3 upper triangular
static
inline void mayo_5_P3_times_S2_neon(const uint64_t *_P3, uint8x16_t *S2_multabs, uint64_t *_acc){
    size_t k,c;
    const uint8x16_t *P3 = (uint8x16_t *) _P3;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );

    size_t cols_used = 0;
    for (size_t r = 0; r < O_MAX; r++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*4] = {0};

        for (c=r; c < O_MAX; c++)
        {
            uint8x16_t in_odd0 = P3[4*cols_used];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = P3[4*cols_used + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = P3[4*cols_used + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = P3[4*cols_used + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;
            cols_used ++;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[8*k]     ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd0);
                temp[8*k + 1] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even0);
                temp[8*k + 2] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd1);
                temp[8*k + 3] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even1);
                temp[8*k + 4] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd2);
                temp[8*k + 5] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even2);
                temp[8*k + 6] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_odd3);
                temp[8*k + 7] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*c + k], in_even3);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k + 1 < K_MAX; k+=2)
        {
            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k] ^= temp[4*k] ^ (t0 << 4);
            acc[4*(r*K_MAX) + 4*k + 4] ^= temp[4*k+1] ^ t0;

            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k+2] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 1] ^= temp[4*k+2] ^ (t1 << 4);
            acc[4*(r*K_MAX) + 4*k + 5] ^= temp[4*k+3] ^ t1;

            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k+4] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 2] ^= temp[4*k+4] ^ (t2 << 4);
            acc[4*(r*K_MAX) + 4*k + 6] ^= temp[4*k+5] ^ t2;

            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k+6] >> 4)) & low_nibble_mask;
            acc[4*(r*K_MAX) + 4*k + 3] ^= temp[4*k+6] ^ (t3 << 4);
            acc[4*(r*K_MAX) + 4*k + 7] ^= temp[4*k+7] ^ t3;
        }
    }
}


static
inline void mayo_5_S1t_times_PS1_neon(const uint64_t *_PS1, uint8x16_t *S1_multabs, uint64_t *_acc){
    mayo_5_Vt_times_Pv_neon(_PS1, S1_multabs, _acc);
}

static
inline void mayo_5_S2t_times_PS2_neon(const uint64_t *_PS2, uint8x16_t *S2_multabs, uint64_t *_acc){
    const uint8x16_t *PS2 = (uint8x16_t *) _PS2;
    uint8x16_t *acc = (uint8x16_t *) _acc;
    const uint8x16_t low_nibble_mask  = vdupq_n_u8( 0xf );
    size_t k;

    for (size_t c = 0; c < K_MAX; c++)
    {
        // do multiplications for one row and accumulate results in temporary format
        uint8x16_t temp[K_OVER_2*2*4] = {0};
        for (size_t r = 0; r < O_MAX; r++)
        {
            uint8x16_t in_odd0 = PS2[4*r*K_MAX + 4*c];
            uint8x16_t in_even0 = (in_odd0 >> 4) & low_nibble_mask;
            in_odd0 &= low_nibble_mask;
            uint8x16_t in_odd1 = PS2[4*r*K_MAX + 4*c + 1];
            uint8x16_t in_even1 = (in_odd1 >> 4) & low_nibble_mask;
            in_odd1 &= low_nibble_mask;
            uint8x16_t in_odd2 = PS2[4*r*K_MAX + 4*c + 2];
            uint8x16_t in_even2 = (in_odd2 >> 4) & low_nibble_mask;
            in_odd2 &= low_nibble_mask;
            uint8x16_t in_odd3 = PS2[4*r*K_MAX + 4*c + 3];
            uint8x16_t in_even3 = (in_odd3 >> 4) & low_nibble_mask;
            in_odd3 &= low_nibble_mask;

            for (k = 0; k < K_OVER_2; k++)
            {
                temp[8*k]     ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_odd0);
                temp[8*k + 1] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_even0);
                temp[8*k + 2] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_odd1);
                temp[8*k + 3] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_even1);
                temp[8*k + 4] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_odd2);
                temp[8*k + 5] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_even2);
                temp[8*k + 6] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_odd3);
                temp[8*k + 7] ^= vqtbl1q_u8(S2_multabs[K_OVER_2*r + k], in_even3);
            }
        }

        // convert to normal format and add to accumulator
        for (k = 0; k+1 < K_MAX; k+=2)
        {
            uint8x16_t t0 = (temp[4*k + 1] ^ (temp[4*k] >> 4)) & low_nibble_mask;
            acc[4*(k*K_MAX) + 4*c] ^= temp[4*k] ^ (t0 << 4);
            acc[4*((k+1)*K_MAX) + 4*c] ^= temp[4*k+1] ^ t0;

            uint8x16_t t1 = (temp[4*k + 3] ^ (temp[4*k+2] >> 4)) & low_nibble_mask;
            acc[4*(k*K_MAX) + 4*c + 1] ^= temp[4*k+2] ^ (t1 << 4);
            acc[4*((k+1)*K_MAX) + 4*c + 1] ^= temp[4*k+3] ^ t1;

            uint8x16_t t2 = (temp[4*k + 5] ^ (temp[4*k+4] >> 4)) & low_nibble_mask;
            acc[4*(k*K_MAX) + 4*c + 2] ^= temp[4*k+4] ^ (t2 << 4);
            acc[4*((k+1)*K_MAX) + 4*c + 2] ^= temp[4*k+5] ^ t2;

            uint8x16_t t3 = (temp[4*k + 7] ^ (temp[4*k+6] >> 4)) & low_nibble_mask;
            acc[4*(k*K_MAX) + 4*c + 3] ^= temp[4*k+6] ^ (t3 << 4);
            acc[4*((k+1)*K_MAX) + 4*c + 3] ^= temp[4*k+7] ^ t3;
        }
    }
}


#undef K_OVER_2
#endif

