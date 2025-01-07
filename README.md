# MAYO-C

![MAYO-C workflow](https://github.com/PQCMayo/MAYO-C/actions/workflows/cmake.yml/badge.svg)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

```
This code is part of a NIST submission for the PQC signatures call.
```

MAYO-C is a C library implementation of [MAYO](https://pqmayo.org), a multivariate quadratic signature scheme. It implements the following parameter sets:

| Parameter Set | NIST Security Level | n | m | o | k | q | sk size | pk size | sig size |
| --- | ---- | -- | -- | -- | -- | -- | -- | -- | -- |
| MAYO_1 | 1 | 86 | 78 | 8 | 10 | 16 | 24 B | 1420 B | 454 B |
| MAYO_2 | 1 | 81 | 64 | 17 | 4 | 16 | 24 B | 4912 B | 186 B |
| MAYO_3 | 3 | 118 | 108 | 10 | 11 | 16 | 32 B | 2986 B | 681 B |
| MAYO_5 | 5 | 154 | 142 | 12 | 12 | 16 | 40 B | 5554 B | 964 B |

## Requirements

- CMake (version 3.10 or later)
- C99-compatible compiler
- Valgrind (for dynamic testing)
- Clang static analyzer (version 10 or later, for static analysis)

## Build

- `mkdir -p build`
- `cd build`
- `cmake <Build Options> ..`
- `make`

The following build options have been used to report performance numbers in the specification:

1. Reference: `cmake -DMAYO_BUILD_TYPE=ref -DENABLE_AESNI=OFF ..`
2. Optimized (AES-NI enabled): `cmake -DMAYO_BUILD_TYPE=opt -DENABLE_AESNI=ON ..`
3. Optimized (AES-NI disabled): `cmake -DMAYO_BUILD_TYPE=opt -DENABLE_AESNI=OFF ..`
4. AVX2: `cmake -DMAYO_BUILD_TYPE=avx2 -DENABLE_AESNI=ON ..`
5. A64 M1/M2/M3 NEON: `cmake -DMAYO_BUILD_TYPE=neon -DENABLE_AESNEON=ON ..`
6. A64 RPi4 Cortex-A72 NEON: `cmake -DMAYO_BUILD_TYPE=neon -DENABLE_AESNEON=OFF ..`

## Build options

CMake build options can be specified with `-D<BUILD_OPTION>=<VALUE>`.

### ENABLE_TESTS

Builds a test harness for the library, the default value is `ON`.

### ENABLE_CT_TESTING

Builds the library with instrumentation for constant-time behavior testing, the default value is `OFF`. Valgrind development files are used for this build option.

### ENABLE_PARAMS_DYNAMIC

Builds the library as a single library dynamically supporting all mayo parameter sets. If the option is turned off, multiple libraries for each parameter sets are built, which usually comes with a performance gain. The default value is `OFF`.

### ENABLE_STRICT

Builds the library in strict mode: warnings terminate compilation). The default value is `ON`.

### MAYO_BUILD_TYPE

Specifies the build type for which Mayo is built. The options are `ref`, `opt` and `avx2`. The effect is the following:

- `ref` builds MAYO implemented with portable C code for native target architecture, using run-time parameters.
- `opt` builds MAYO implemented with optimized portable C code, compiled with `-march=native` (if available) and AES acceleration (if available)
- `avx2` builds MAYO implemented with optimized AVX2 code, compiled with `-march=native` (if available) and AES acceleration (if available)
- `neon` builds MAYO implemented with optimized NEON code, compiled with `-march=native` (if available) and AES acceleration (if available)

The default build type if none is specified is `opt`.

### CMAKE_BUILD_TYPE

Can be used to specify special build types. The options are:

- `ASAN`: Builds with AddressSanitizer memory error detector.
- `MSAN`: Builds with MemorySanitizer detector for uninitialized reads.
- `LSAN`: Builds with LeakSanitizer for run-time memory leak detection.
- `UBSAN`: Builds with UndefinedBehaviorSanitizer for undefined behavior detection.

The default build type uses `-O3 -Wstrict-prototypes -Wno-error=strict-prototypes -fvisibility=hidden -Wno-error=implicit-function-declaration -Wno-error=attributes`.

## Build artifacts

The following artifacts are built:

- `libmayo_common_sys.a`: library with common crypto - AES, Keccak and system random number generator.
- `libmary_common_test.a`: library with common crypto for deterministic tests - AES, Keccak and CTR-DRBG PRNG.
- `libmayo_<level>.a`: library for `MAYO_<level>`.
- `libmayo_<level>_test`: library for `MAYO_<level>`, only for test, using the deterministic CTR-DRBG as backend.
- `libmayo_<level>_nistapi.a`: library for `MAYO_<level>` against the NIST API.
- `libmayo_<level>_nistapi_test.a`: library for `MAYO_<level>` against the NIST API. Only for test, using the deterministic CTR-DRBG as backend.
- `mayo_bench_<param>`: Benchmarking suites.
- `mayo_test_kat_<param>` (`opt`, `avx2` and `neon`), `mayo_test_kat` (`ref`): KAT test suites.
- `mayo_test_scheme_<param>` (`opt`, `avx2` and `neon`), `mayo_test_scheme` (`ref`): Self-test suites.
- `PQCgenKAT_sign_<param>`: App for generating NIST KAT.
- `example_<param>` (`opt`, `avx2` and `neon`), `example_mayo` (`ref`): Example app using the MAYO API.
- `example_nistapi_<param>`: Example app using the NIST API.

## Test

In the `build` directory, run: `make test`.

The test harness consists of the following units:

- KAT test: tests against the KAT files in the `KAT` folder - `MAYO_<level>_KAT`
- Self-tests: runs random self-tests (key-generation, signing and verifying) - `MAYO_<level>_SELFTEST`

## Known Answer Tests (KAT)

KAT are available in folder `KAT`. They can be generated by running the apps built in the `apps` folder:

- `apps/PQCgenKAT_sign_mayo_1`
- `apps/PQCgenKAT_sign_mayo_2`
- `apps/PQCgenKAT_sign_mayo_3`
- `apps/PQCgenKAT_sign_mayo_5`

A successful execution will generate the `.req` and `.rsp` files.

KAT verification is done as part of the test harness (see previous section).

## Benchmarks

A benchmarking suite is built and runs with the following command, where `params` specifies the MAYO parameter set and `runs` the number of benchmark runs:

If `MAYO_BUILD_TYPE` is `opt`, `avx2` or `neon`:
- `test/mayo_bench_<param> <runs>`,
- On Apple M1/M2/M3 chips, this must be run with root (`sudo`) permission.

If `MAYO_BUILD_TYPE` is `ref`:
- `test/mayo_bench <param> <runs>`,

The benchmarks profile the `MAYO.CompactKeyGen`, `MAYO.expandSK`, `MAYO.expandSK`, `MAYO.sign` and `MAYO.verify` functions. The results are reported in CPU cycles if available on the host platform, and timing in nanoseconds otherwise.

## Examples

Example code that demonstrates how to use MAYO both via the MAYO API and NIST API are available in the `apps` folder:

- `apps/example.c`: Example with the MAYO API.
- `apps/example_nistapi.c`: Example with the NIST API.

## Project Structure

- `apps`: Applications: KAT generation application
- `include`: MAYO public header files
- `KAT`: Known Answer Test files
- `src`: MAYO source code
- `src/mayo_<x>` MAYO implementation with NIST signature API
- `src/common`: MAYO common components (RNG, AES, SHAKE)
- `src/generic`: MAYO generic C source code
- `src/<arch>`: MAYO <arch> specific source code
- `test`: MAYO test code

## License

MAYO-C is licensed under Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

Third party code is used in some test and common code files:

- `common/aes_c.c`; MIT: "Copyright (c) 2016 Thomas Pornin <pornin@bolet.org>"
- `common/aes128ctr.c`: MIT: "Copyright (c) 2016-2021 Open Quantum Safe project" and Public Domain
- `common/aes_neon.c`: MIT: "Copyright (c) 2024 ChristerKnorborg" and Public Domain
- `common/fips202.c`: Public Domain
- `common/randombytes_system.c`: MIT: Copyright (c) 2017 Daan Sprenkels <hello@dsprenkels.com>
- `apps/PQCgenKAT_sign.c`, `common/randombytes_ctrdrbg.c`, `test/test_kat.c`: by NIST (Public Domain)
- `test/m1cycles.{c,h}`, Apache 2.0 and Public Domain

See also the SPDX License Identifiers in the respective files.

## Citing

Bibtext:

```
@manual{mayo-c,
  title        = {MAYO C implementation},
  author       = {Ward Beullens and Fabio Campos and Sof\'{i}a Celi and Basil Hess and Matthias J. Kannwischer},
  note         = {Available at \url{https://github.com/PQCMayo/MAYO-C}. Accessed June, 2023},
  month        = jun,
  year         = {2023}
}
```
