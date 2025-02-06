# SPDX-License-Identifier: Apache-2.0

# Script go generate constants for the parameter sets.

import math

DEFAULT_PARAMETERS = {
    "MAYO_1": {
        "name": "mayo1",
        "n": 86,
        "m": 78,
        "o": 8,
        "k": 10,
        "q": 16,
        "pk_seed_bytes": 16,
        "sk_seed_bytes": 24,
        "salt_bytes": 24,
        "digest_bytes": 32
    },
    "MAYO_2": {
        "name": "mayo2",
        "n": 81,
        "m": 64,
        "o": 17,
        "k": 4,
        "q": 16,
        "pk_seed_bytes": 16,
        "sk_seed_bytes": 24,
        "salt_bytes": 24,
        "digest_bytes": 32
    },
    "MAYO_3": {
        "name": "mayo3",
        "n": 118,
        "m": 108,
        "o": 10,
        "k": 11,
        "q": 16,
        "pk_seed_bytes": 16,
        "sk_seed_bytes": 32,
        "salt_bytes": 32,
        "digest_bytes": 48
    },
    "MAYO_5": {
        "name": "mayo5",
        "n": 154,
        "m": 142,
        "o": 12,
        "k": 12,
        "q": 16,
        "pk_seed_bytes": 16,
        "sk_seed_bytes": 40,
        "salt_bytes": 40,
        "digest_bytes": 64
    },
}

for param in DEFAULT_PARAMETERS:
    n = DEFAULT_PARAMETERS[param]["n"]
    m = DEFAULT_PARAMETERS[param]["m"]
    o = DEFAULT_PARAMETERS[param]["o"]
    k = DEFAULT_PARAMETERS[param]["k"]
    q = DEFAULT_PARAMETERS[param]["q"]
    pk_seed_bytes = DEFAULT_PARAMETERS[param]["pk_seed_bytes"]
    sk_seed_bytes = DEFAULT_PARAMETERS[param]["sk_seed_bytes"]
    salt_bytes = DEFAULT_PARAMETERS[param]["salt_bytes"]
    digest_bytes = DEFAULT_PARAMETERS[param]["digest_bytes"]

    v = n - o
    q_bytes = (math.log(q, 2)/8)
    m_bytes = math.ceil(q_bytes*m)
    O_bytes = math.ceil((n - o)*o*q_bytes)
    v_bytes = math.ceil((n - o)*q_bytes)
    r_bytes = math.ceil(k*o*q_bytes)
    P1_bytes = math.ceil(m*math.comb((n-o+1), 2)*q_bytes)
    P2_bytes = math.ceil(m*(n - o)*o*q_bytes)
    P3_bytes = math.ceil(m*math.comb((o+1), 2)*q_bytes)
    m_vec_limbs = math.ceil(m/16)

    R_bytes = salt_bytes
    sig_bytes = math.ceil(k * n * q_bytes) + salt_bytes
    epk_bytes = P1_bytes + P2_bytes + P3_bytes
    cpk_bytes = P3_bytes + pk_seed_bytes
    csk_bytes = sk_seed_bytes
    esk_bytes = sk_seed_bytes + O_bytes + P1_bytes + P2_bytes

    DEFAULT_PARAMETERS[param]["v"] = v
    DEFAULT_PARAMETERS[param]["O_bytes"] = O_bytes
    DEFAULT_PARAMETERS[param]["v_bytes"] = v_bytes
    DEFAULT_PARAMETERS[param]["r_bytes"] = r_bytes
    DEFAULT_PARAMETERS[param]["P1_bytes"] = P1_bytes
    DEFAULT_PARAMETERS[param]["P2_bytes"] = P2_bytes
    DEFAULT_PARAMETERS[param]["P3_bytes"] = P3_bytes
    DEFAULT_PARAMETERS[param]["sig_bytes"] = sig_bytes
    DEFAULT_PARAMETERS[param]["cpk_bytes"] = cpk_bytes
    DEFAULT_PARAMETERS[param]["csk_bytes"] = csk_bytes
    DEFAULT_PARAMETERS[param]["m_bytes"] = m_bytes
    DEFAULT_PARAMETERS[param]["pk_seed_bytes"] = pk_seed_bytes
    DEFAULT_PARAMETERS[param]["sk_seed_bytes"] = sk_seed_bytes
    DEFAULT_PARAMETERS[param]["salt_bytes"] = salt_bytes
    DEFAULT_PARAMETERS[param]["digest_bytes"] = digest_bytes
    DEFAULT_PARAMETERS[param]["m_vec_limbs"] = m_vec_limbs

    print("#define " + param + "_n " + str(n))
    print("#define " + param + "_m " + str(m))
    print("#define " + param + "_m_vec_limbs " + str(m_vec_limbs))
    print("#define " + param + "_o " + str(o))
    print("#define " + param + "_v " + str(v))
    print("#define " + param + "_A_cols (" + param + "_k * " + param + "_o + 1)")
    print("#define " + param + "_k " + str(k))
    print("#define " + param + "_q " + str(q))
    print("#define " + param + "_m_bytes " + str(m_bytes))
    print("#define " + param + "_O_bytes " + str(O_bytes))
    print("#define " + param + "_v_bytes " + str(v_bytes))
    print("#define " + param + "_r_bytes " + str(r_bytes))
    print("#define " + param + "_P1_bytes " + str(P1_bytes))
    print("#define " + param + "_P2_bytes " + str(P2_bytes))
    print("#define " + param + "_P3_bytes " + str(P3_bytes))
    print("#define " + param + "_csk_bytes " + str(csk_bytes))
    print("#define " + param + "_cpk_bytes " + str(cpk_bytes))
    print("#define " + param + "_sig_bytes " + str(sig_bytes))
    print("#define " + param + "_f_tail F_TAIL_" + str(m))
    print("#define " + param + "_f_tail_arr f_tail_" + str(m))
    print("#define " + param + "_salt_bytes " + str(salt_bytes))
    print("#define " + param + "_digest_bytes " + str(digest_bytes))
    print("#define " + param + "_pk_seed_bytes " + str(pk_seed_bytes))
    print("#define " + param + "_sk_seed_bytes " + str(sk_seed_bytes))



    print("\n")

maxvals = {}

for param in DEFAULT_PARAMETERS:
    for paramval in DEFAULT_PARAMETERS[param]:
        if paramval not in maxvals or maxvals[paramval] < DEFAULT_PARAMETERS[param][paramval]:
            maxvals[paramval] = DEFAULT_PARAMETERS[param][paramval]

for max in maxvals:
    print("#define " + max.upper() + "_MAX " + str(maxvals[max]))