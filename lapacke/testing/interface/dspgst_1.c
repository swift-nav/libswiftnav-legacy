/*****************************************************************************
  Copyright (c) 2010, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/
/*  Contents: test routine for C interface to LAPACK
*   Author: Intel Corporation
*   Created in March, 2010
*
* Purpose
*
* dspgst_1 is the test program for the C interface to LAPACK
* routine dspgst
* The program doesn't require an input, the input data is hardcoded in the
* test program.
* The program tests the C interface in the four combinations:
*   1) column-major layout, middle-level interface
*   2) column-major layout, high-level interface
*   3) row-major layout, middle-level interface
*   4) row-major layout, high-level interface
* The output of the C interface function is compared to those obtained from
* the corresponiding LAPACK routine with the same input data, and the
* comparison diagnostics is then printed on the standard output having PASSED
* keyword if the test is passed, and FAILED keyword if the test isn't passed.
*****************************************************************************/
#include <stdio.h>
#include "lapacke.h"
#include "lapacke_utils.h"
#include "test_utils.h"

static void init_scalars_dspgst( lapack_int *itype, char *uplo, lapack_int *n );
static void init_ap( lapack_int size, double *ap );
static void init_bp( lapack_int size, double *bp );
static int compare_dspgst( double *ap, double *ap_i, lapack_int info,
                           lapack_int info_i, lapack_int n );

int main(void)
{
    /* Local scalars */
    lapack_int itype, itype_i;
    char uplo, uplo_i;
    lapack_int n, n_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *ap = NULL, *ap_i = NULL;
    double *bp = NULL, *bp_i = NULL;
    double *ap_save = NULL;
    double *ap_r = NULL;
    double *bp_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dspgst( &itype, &uplo, &n );
    itype_i = itype;
    uplo_i = uplo;
    n_i = n;

    /* Allocate memory for the LAPACK routine arrays */
    ap = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );
    bp = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    ap_i = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );
    bp_i = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );

    /* Allocate memory for the backup arrays */
    ap_save = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    ap_r = (double *)LAPACKE_malloc( n*(n+1)/2 * sizeof(double) );
    bp_r = (double *)LAPACKE_malloc( n*(n+1)/2 * sizeof(double) );

    /* Initialize input arrays */
    init_ap( (n*(n+1)/2), ap );
    init_bp( (n*(n+1)/2), bp );

    /* Backup the ouptut arrays */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_save[i] = ap[i];
    }

    /* Call the LAPACK routine */
    dspgst_( &itype, &uplo, &n, ap, bp, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        bp_i[i] = bp[i];
    }
    info_i = LAPACKE_dspgst_work( LAPACK_COL_MAJOR, itype_i, uplo_i, n_i, ap_i,
                                  bp_i );

    failed = compare_dspgst( ap, ap_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dspgst\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dspgst\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        bp_i[i] = bp[i];
    }
    info_i = LAPACKE_dspgst( LAPACK_COL_MAJOR, itype_i, uplo_i, n_i, ap_i,
                             bp_i );

    failed = compare_dspgst( ap, ap_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dspgst\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dspgst\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        bp_i[i] = bp[i];
    }

    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, n, bp_i, bp_r );
    info_i = LAPACKE_dspgst_work( LAPACK_ROW_MAJOR, itype_i, uplo_i, n_i, ap_r,
                                  bp_r );

    LAPACKE_dpp_trans( LAPACK_ROW_MAJOR, uplo, n, ap_r, ap_i );

    failed = compare_dspgst( ap, ap_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dspgst\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dspgst\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        bp_i[i] = bp[i];
    }

    /* Init row_major arrays */
    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, n, bp_i, bp_r );
    info_i = LAPACKE_dspgst( LAPACK_ROW_MAJOR, itype_i, uplo_i, n_i, ap_r,
                             bp_r );

    LAPACKE_dpp_trans( LAPACK_ROW_MAJOR, uplo, n, ap_r, ap_i );

    failed = compare_dspgst( ap, ap_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dspgst\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dspgst\n" );
    }

    /* Release memory */
    if( ap != NULL ) {
        LAPACKE_free( ap );
    }
    if( ap_i != NULL ) {
        LAPACKE_free( ap_i );
    }
    if( ap_r != NULL ) {
        LAPACKE_free( ap_r );
    }
    if( ap_save != NULL ) {
        LAPACKE_free( ap_save );
    }
    if( bp != NULL ) {
        LAPACKE_free( bp );
    }
    if( bp_i != NULL ) {
        LAPACKE_free( bp_i );
    }
    if( bp_r != NULL ) {
        LAPACKE_free( bp_r );
    }

    return 0;
}

/* Auxiliary function: dspgst scalar parameters initialization */
static void init_scalars_dspgst( lapack_int *itype, char *uplo, lapack_int *n )
{
    *itype = 1;
    *uplo = 'L';
    *n = 4;

    return;
}

/* Auxiliary functions: dspgst array parameters initialization */
static void init_ap( lapack_int size, double *ap ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ap[i] = 0;
    }
    ap[0] = 2.39999999999999990e-001;
    ap[1] = 3.90000000000000010e-001;
    ap[2] = 4.19999999999999980e-001;
    ap[3] = -1.60000000000000000e-001;
    ap[4] = -1.10000000000000000e-001;
    ap[5] = 7.90000000000000040e-001;
    ap[6] = 6.30000000000000000e-001;
    ap[7] = -2.50000000000000000e-001;
    ap[8] = 4.79999999999999980e-001;
    ap[9] = -2.99999999999999990e-002;
}
static void init_bp( lapack_int size, double *bp ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        bp[i] = 0;
    }
    bp[0] = 2.03960780543711410e+000;
    bp[1] = -1.52970585407783540e+000;
    bp[2] = 2.74562589193457660e-001;
    bp[3] = -4.90290337845460110e-002;
    bp[4] = 1.64012194668567270e+000;
    bp[5] = -2.49981411948373810e-001;
    bp[6] = 6.18856422262437760e-001;
    bp[7] = 7.88748805574805310e-001;
    bp[8] = 6.44266130231023420e-001;
    bp[9] = 6.16063337578069970e-001;
}

/* Auxiliary function: C interface to dspgst results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dspgst( double *ap, double *ap_i, lapack_int info,
                           lapack_int info_i, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        failed += compare_doubles(ap[i],ap_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
