// SPDX-License-Identifier: Apache-2.0

#ifndef M1CYCLES_H
#define M1CYCLES_H

#ifdef TARGET_ARM64

void setup_rdtsc(void);
unsigned long long int rdtsc(void);

#endif

#endif
