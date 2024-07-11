// SPDX-License-Identifier: Apache-2.0 and Unknown

/*
NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <randombytes.h>
#include <mayo.h>

#define MAX_MARKER_LEN         50

#define KAT_SUCCESS             0
#define KAT_FILE_OPEN_ERROR    -1
#define KAT_DATA_ERROR         -3
#define KAT_CRYPTO_FAILURE     -4
#define KAT_VERIFICATION_ERROR -5

static int      FindMarker(FILE *infile, const char *marker);
static int      ReadHex(FILE *infile, unsigned char *A, int Length, char *str);
static int      test_sig_kat(const mayo_params_t *p);

int main(int argc, char *argv[]) {
    int rc = 0;

#ifdef ENABLE_PARAMS_DYNAMIC
    if (argc < 2) {
        printf("One argument needed\n");
        rc = 1;
        goto end;
    }

    if (!strcmp(argv[1], "MAYO-1")) {
        rc = test_sig_kat(&MAYO_1);
    } else if (!strcmp(argv[1], "MAYO-2")) {
        rc = test_sig_kat(&MAYO_2);
    } else if (!strcmp(argv[1], "MAYO-3")) {
        rc = test_sig_kat(&MAYO_3);
    } else if (!strcmp(argv[1], "MAYO-5")) {
        rc = test_sig_kat(&MAYO_5);
    } else {
        return -1;
    }
end:
#else
    rc = test_sig_kat(0);
#endif
    return rc;
}


static int test_sig_kat(const mayo_params_t *p) {
    unsigned char       seed[48];
    unsigned char       *m, *sm, *m1, *sm_rsp;
    size_t  mlen, smlen, mlen1;
    int                 count;
    int                 done;
    unsigned char       pk[PARAM_cpk_bytes(p)], sk[PARAM_csk_bytes(p)];
    int                 ret_val;

    char                fn_rsp[64];
    FILE                *fp_rsp;
    unsigned char       pk_rsp[PARAM_cpk_bytes(p)], sk_rsp[PARAM_csk_bytes(p)];

    sprintf(fn_rsp, "../../KAT/PQCsignKAT_%d_%s.rsp", PARAM_csk_bytes(p), PARAM_name(p));
    if ( (fp_rsp = fopen(fn_rsp, "r")) == NULL ) {
        printf("Couldn't open <%s> for read\n", fn_rsp);
        return KAT_FILE_OPEN_ERROR;
    }

    done = 0;
    do {
        if ( FindMarker(fp_rsp, "count = ") ) {
            ret_val = fscanf(fp_rsp, "%d", &count);
        } else {
            done = 1;
            break;
        }

        if ( !ReadHex(fp_rsp, seed, 48, "seed = ") ) {
            printf("ERROR: unable to read 'seed' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }

        randombytes_init(seed, NULL, 256);

        if ( FindMarker(fp_rsp, "mlen = ") ) {
            ret_val = fscanf(fp_rsp, "%zu", &mlen);
        } else {
            printf("ERROR: unable to read 'mlen' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }

        m = (unsigned char *)calloc(mlen, sizeof(unsigned char));
        m1 = (unsigned char *)calloc(mlen, sizeof(unsigned char));
        sm = (unsigned char *)calloc(mlen + PARAM_sig_bytes(p), sizeof(unsigned char));
        sm_rsp = (unsigned char *)calloc(mlen + PARAM_sig_bytes(p), sizeof(unsigned char));

        if ( !ReadHex(fp_rsp, m, (int)mlen, "msg = ") ) {
            printf("ERROR: unable to read 'msg' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }

        // Generate the public/private keypair
        if ( (ret_val = mayo_keypair(p, pk, sk)) != 0) {
            printf("crypto_sign_keypair returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        if ( !ReadHex(fp_rsp, pk_rsp, PARAM_cpk_bytes(p), "pk = ") ) {
            printf("ERROR: unable to read 'pk' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }
        if ( !ReadHex(fp_rsp, sk_rsp, PARAM_csk_bytes(p), "sk = ") ) {
            printf("ERROR: unable to read 'sk' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }

        if (memcmp(pk, pk_rsp, PARAM_cpk_bytes(p)) != 0) {
            printf("ERROR: pk is different from <%s>\n", fn_rsp);
            return KAT_VERIFICATION_ERROR;
        }
        if (memcmp(sk, sk_rsp, PARAM_csk_bytes(p)) != 0) {
            printf("ERROR: sk is different from <%s>\n", fn_rsp);
            return KAT_VERIFICATION_ERROR;
        }

        if ( (ret_val = mayo_sign(p, sm, &smlen, m, mlen, sk)) != 0) {
            printf("crypto_sign returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }

        if ( !ReadHex(fp_rsp, sm_rsp, smlen, "sm = ") ) {
            printf("ERROR: unable to read 'sm' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }

        if (memcmp(sm, sm_rsp, smlen) != 0) {
            printf("ERROR: sm is different from <%s>\n", fn_rsp);
            return KAT_VERIFICATION_ERROR;
        }


        if ( (ret_val = mayo_open(p, m1, &mlen1, sm, smlen, pk)) != 0) {
            printf("crypto_sign_open returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }

        if ( mlen != mlen1 ) {
            printf("crypto_sign_open returned bad 'mlen': Got <%zu>, expected <%zu>\n", mlen1, mlen);
            return KAT_CRYPTO_FAILURE;
        }

        if ( memcmp(m, m1, mlen) ) {
            printf("crypto_sign_open returned bad 'm' value\n");
            return KAT_CRYPTO_FAILURE;
        }

        free(m);
        free(m1);
        free(sm);
        free(sm_rsp);

    } while ( !done );

    fclose(fp_rsp);

    printf("Known Answer Tests PASSED. \n");
    printf("\n\n");

    return KAT_SUCCESS;
}


//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
static int
FindMarker(FILE *infile, const char *marker) {
    char    line[MAX_MARKER_LEN];
    int     i, len;
    int curr_line;

    len = (int)strlen(marker);
    if ( len > MAX_MARKER_LEN - 1 ) {
        len = MAX_MARKER_LEN - 1;
    }

    for ( i = 0; i < len; i++ ) {
        curr_line = fgetc(infile);
        line[i] = curr_line;
        if (curr_line == EOF ) {
            return 0;
        }
    }
    line[len] = '\0';

    while ( 1 ) {
        if ( !strncmp(line, marker, len) ) {
            return 1;
        }

        for ( i = 0; i < len - 1; i++ ) {
            line[i] = line[i + 1];
        }
        curr_line = fgetc(infile);
        line[len - 1] = curr_line;
        if (curr_line == EOF ) {
            return 0;
        }
        line[len] = '\0';
    }

    // shouldn't get here
    return 0;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
static int
ReadHex(FILE *infile, unsigned char *A, int Length, char *str) {
    int         i, ch, started;
    unsigned char   ich;

    if ( Length == 0 ) {
        A[0] = 0x00;
        return 1;
    }
    memset(A, 0x00, Length);
    started = 0;
    if ( FindMarker(infile, str) )
        while ( (ch = fgetc(infile)) != EOF ) {
            if ( !isxdigit(ch) ) {
                if ( !started ) {
                    if ( ch == '\n' ) {
                        break;
                    } else {
                        continue;
                    }
                } else {
                    break;
                }
            }
            started = 1;
            if ( (ch >= '0') && (ch <= '9') ) {
                ich = ch - '0';
            } else if ( (ch >= 'A') && (ch <= 'F') ) {
                ich = ch - 'A' + 10;
            } else if ( (ch >= 'a') && (ch <= 'f') ) {
                ich = ch - 'a' + 10;
            } else { // shouldn't ever get here
                ich = 0;
            }

            for ( i = 0; i < Length - 1; i++ ) {
                A[i] = (A[i] << 4) | (A[i + 1] >> 4);
            }
            A[Length - 1] = (A[Length - 1] << 4) | ich;
        } else {
        return 0;
    }

    return 1;
}

