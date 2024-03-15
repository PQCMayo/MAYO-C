# Nibbling-MAYO Artifact

AVX2 implementations of MAYO, a multivariate quadratic signature scheme, as described in the paper **Nibbling MAYO: Optimized Implementations for AVX2 and Cortex-M4** available [here](https://eprint.iacr.org/2023/1683.pdf).

It consists of two variants:

1. A slower version compatible with the round-1 specification of MAYO using bitsliced representation. This version is contained in the [main](https://github.com/PQCMayo/MAYO-C/tree/main) branch of this repository.
2. A faster version that changes representation of keys and PRNG output to nibble-sliced representation. This version is contained in the [nibbling-mayo branch](https://github.com/PQCMayo/MAYO-C/tree/nibbling-mayo) of this repository.

All implementations implement the following parameter sets:

| Parameter Set | NIST Security Level | n | m | o | k | q | sk size | pk size | sig size |
| --- | ---- | -- | -- | -- | -- | -- | -- | -- | -- |
| MAYO_1 | 1 | 66 | 64 | 8 | 9 | 16 | 24 B | 1168 B | 321 B |
| MAYO_2 | 1 | 78 | 64 | 18 | 4 | 16 | 24 B | 5488 B | 180 B |
| MAYO_3 | 3 | 99 | 96 | 10 | 11 | 16 | 32 B | 2656 B | 577 B |
| MAYO_5 | 5 | 133 | 128 | 12 | 12 | 16 | 40 B | 5008 B | 838 B |

This file contains information about the build environment, build instructions and benchmarking programs to reproduce the results of the paper.

## Build Requirements

- CMake (version 3.5 or later)
- C99 compatible compiler
- Intel Haswell or later for AVX2 optimizations

## Evaluation environment used in the paper

The following CPUs were used:
- Intel Xeon X3-1245 v5 (Skylake)
- Intel Xeon Gold 6338 (Ice Lake)

Operating system: Ubuntu 22.04.3 LTS

Compiler: Ubuntu clang version 14.0.0-1ubuntu1.1

Turbo Boost was disabled in the UEFI.

## Build

The following steps build MAYO-C with AVX2 optimizations using clang:

- `git checkout main` or `git checkout nibbling-mayo`
- `mkdir -p build`
- `cd build`
- `cmake -DMAYO_BUILD_TYPE=avx2 -DCMAKE_C_COMPILER=clang ..`
- `make -j`

Other options for `MAYO_BUILD_TYPE` are `opt` and `ref`.

## Running benchmarks

Benchmarking programs are available in folder `test` after a successful build. To reproduce the results of Table 1, run:

- `test/mayo_bench_MAYO_1 <repetitions>`
- `test/mayo_bench_MAYO_2 <repetitions>`
- `test/mayo_bench_MAYO_3 <repetitions>`
- `test/mayo_bench_MAYO_5 <repetitions>`

The benchmarks report the median and the average cycle count from the specified number of repetitions (10000 in the paper), for the following MAYO operations:

- KeyGen (mayo_keygen)
- ExpandSK (mayo_expand_sk)
- ExpandPK (mayo_expand_pk)
- ExpandSK + Sign (mayo_sign)
- ExpandPK + Verify (mayo_verify)

To reproduce the results of Table 3, run:

- `test/mayo_bench_table3_MAYO_1 <repetitions>`
- `test/mayo_bench_table3_MAYO_2 <repetitions>`
- `test/mayo_bench_table3_MAYO_3 <repetitions>`
- `test/mayo_bench_table3_MAYO_5 <repetitions>`

## KAT and self-tests

KAT and selftests are available after a successful build. To run them, run one of the following in the build folder:

- `make test` or `ctest`

## Detailed user guide

A more detailed [README](README.md) with all configuration options of the MAYO-C library is available.

## License

The implementations are licenced under Apache-2.0, see the LICENSE and NOTICE files.