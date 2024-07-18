# SPDX-License-Identifier: Apache-2.0

# AddressSanitizer
set(CMAKE_C_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
    FORCE)

# LeakSanitizer
set(CMAKE_C_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during LeakSanitizer builds."
    FORCE)

# MemorySanitizer
set(CMAKE_C_FLAGS_MSAN
    "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during MemorySanitizer builds."
    FORCE)

# UndefinedBehaviour
set(CMAKE_C_FLAGS_UBSAN
    "-fsanitize=undefined"
    CACHE STRING "Flags used by the C compiler during UndefinedBehaviourSanitizer builds."
    FORCE)

set(CMAKE_C_FLAGS_COVERAGE
    "-fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C compiler during Coverage builds."
    FORCE)

# CT-testing configs
set(CMAKE_C_FLAGS_CTOS
    "-Os -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTO0
    "-O0 -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTO1
    "-O1 -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTO2
    "-O2 -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTO3
    "-O3 -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTOSNOVEC
    "-Os -fno-vectorize -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTO0NOVEC
    "-O0 -fno-vectorize -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTO1NOVEC
    "-O1 -fno-vectorize -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTO2NOVEC
    "-O2 -fno-vectorize -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)

set(CMAKE_C_FLAGS_CTO3NOVEC
    "-O3 -fno-vectorize -gdwarf-4"
    CACHE STRING "Flags used by the C compiler during CT builds."
    FORCE)
