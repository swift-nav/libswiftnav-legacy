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
* ztbtrs_1 is the test program for the C interface to LAPACK
* routine ztbtrs
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

static void init_scalars_ztbtrs( char *uplo, char *trans, char *diag,
                                 lapack_int *n, lapack_int *kd,
                                 lapack_int *nrhs, lapack_int *ldab,
                                 lapack_int *ldb );
static void init_ab( lapack_int size, lapack_complex_double *ab );
static void init_b( lapack_int size, lapack_complex_double *b );
static int compare_ztbtrs( lapack_complex_double *b, lapack_complex_double *b_i,
                           lapack_int info, lapack_int info_i, lapack_int ldb,
                           lapack_int nrhs );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
    char trans, trans_i;
    char diag, diag_i;
    lapack_int n, n_i;
    lapack_int kd, kd_i;
    lapack_int nrhs, nrhs_i;
    lapack_int ldab, ldab_i;
    lapack_int ldab_r;
    lapack_int ldb, ldb_i;
    lapack_int ldb_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_double *ab = NULL, *ab_i = NULL;
    lapack_complex_double *b = NULL, *b_i = NULL;
    lapack_complex_double *b_save = NULL;
    lapack_complex_double *ab_r = NULL;
    lapack_complex_double *b_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ztbtrs( &uplo, &trans, &diag, &n, &kd, &nrhs, &ldab, &ldb );
    ldab_r = n+2;
    ldb_r = nrhs+2;
    uplo_i = uplo;
    trans_i = trans;
    diag_i = diag;
    n_i = n;
    kd_i = kd;
    nrhs_i = nrhs;
    ldab_i = ldab;
    ldb_i = ldb;

    /* Allocate memory for the LAPACK routine arrays */
    ab = (lapack_complex_double *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_double) );
    b = (lapack_complex_double *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_double) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (lapack_complex_double *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_double) );
    b_i = (lapack_complex_double *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_double) );

    /* Allocate memory for the backup arrays */
    b_save = (lapack_complex_double *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_double) );

    /* Allocate memory for the row-major arrays */
    ab_r = (lapack_complex_double *)
        LAPACKE_malloc( (kd+1)*(n+2) * sizeof(lapack_complex_double) );
    b_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(nrhs+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_b( ldb*nrhs, b );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_save[i] = b[i];
    }

    /* Call the LAPACK routine */
    ztbtrs_( &uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b_save[i];
    }
    info_i = LAPACKE_ztbtrs_work( LAPACK_COL_MAJOR, uplo_i, trans_i, diag_i,
                                  n_i, kd_i, nrhs_i, ab_i, ldab_i, b_i, ldb_i );

    failed = compare_ztbtrs( b, b_i, info, info_i, ldb, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ztbtrs\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ztbtrs\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b_save[i];
    }
    info_i = LAPACKE_ztbtrs( LAPACK_COL_MAJOR, uplo_i, trans_i, diag_i, n_i,
                             kd_i, nrhs_i, ab_i, ldab_i, b_i, ldb_i );

    failed = compare_ztbtrs( b, b_i, info, info_i, ldb, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ztbtrs\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ztbtrs\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b_save[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, kd+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    info_i = LAPACKE_ztbtrs_work( LAPACK_ROW_MAJOR, uplo_i, trans_i, diag_i,
                                  n_i, kd_i, nrhs_i, ab_r, ldab_r, b_r, ldb_r );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, nrhs, b_r, nrhs+2, b_i, ldb );

    failed = compare_ztbtrs( b, b_i, info, info_i, ldb, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ztbtrs\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ztbtrs\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, kd+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    info_i = LAPACKE_ztbtrs( LAPACK_ROW_MAJOR, uplo_i, trans_i, diag_i, n_i,
                             kd_i, nrhs_i, ab_r, ldab_r, b_r, ldb_r );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, nrhs, b_r, nrhs+2, b_i, ldb );

    failed = compare_ztbtrs( b, b_i, info, info_i, ldb, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ztbtrs\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ztbtrs\n" );
    }

    /* Release memory */
    if( ab != NULL ) {
        LAPACKE_free( ab );
    }
    if( ab_i != NULL ) {
        LAPACKE_free( ab_i );
    }
    if( ab_r != NULL ) {
        LAPACKE_free( ab_r );
    }
    if( b != NULL ) {
        LAPACKE_free( b );
    }
    if( b_i != NULL ) {
        LAPACKE_free( b_i );
    }
    if( b_r != NULL ) {
        LAPACKE_free( b_r );
    }
    if( b_save != NULL ) {
        LAPACKE_free( b_save );
    }

    return 0;
}

/* Auxiliary function: ztbtrs scalar parameters initialization */
static void init_scalars_ztbtrs( char *uplo, char *trans, char *diag,
                                 lapack_int *n, lapack_int *kd,
                                 lapack_int *nrhs, lapack_int *ldab,
                                 lapack_int *ldb )
{
    *uplo = 'L';
    *trans = 'N';
    *diag = 'N';
    *n = 4;
    *kd = 2;
    *nrhs = 2;
    *ldab = 9;
    *ldb = 8;

    return;
}

/* Auxiliary functions: ztbtrs array parameters initialization */
static void init_ab( lapack_int size, lapack_complex_double *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    ab[0] = lapack_make_complex_double( -1.93999999999999990e+000,
                                        4.42999999999999970e+000 );
    ab[9] = lapack_make_complex_double( 4.12000000000000010e+000,
                                        -4.26999999999999960e+000 );
    ab[18] = lapack_make_complex_double( 4.29999999999999990e-001,
                                         -2.66000000000000010e+000 );
    ab[27] = lapack_make_complex_double( 4.40000000000000000e-001,
                                         1.00000000000000010e-001 );
    ab[1] = lapack_make_complex_double( -3.39000000000000010e+000,
                                        3.43999999999999990e+000 );
    ab[10] = lapack_make_complex_double( -1.84000000000000010e+000,
                                         5.53000000000000020e+000 );
    ab[19] = lapack_make_complex_double( 1.74000000000000000e+000,
                                         -4.00000000000000010e-002 );
    ab[28] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    ab[2] = lapack_make_complex_double( 1.62000000000000010e+000,
                                        3.68000000000000020e+000 );
    ab[11] = lapack_make_complex_double( -2.77000000000000000e+000,
                                         -1.92999999999999990e+000 );
    ab[20] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    ab[29] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
}
static void init_b( lapack_int size, lapack_complex_double *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    b[0] = lapack_make_complex_double( -8.85999999999999940e+000,
                                       -3.87999999999999990e+000 );
    b[8] = lapack_make_complex_double( -2.40900000000000000e+001,
                                       -5.26999999999999960e+000 );
    b[1] = lapack_make_complex_double( -1.55700000000000000e+001,
                                       -2.34100000000000000e+001 );
    b[9] = lapack_make_complex_double( -5.79699999999999990e+001,
                                       8.14000000000000060e+000 );
    b[2] = lapack_make_complex_double( -7.62999999999999990e+000,
                                       2.27800000000000010e+001 );
    b[10] = lapack_make_complex_double( 1.90900000000000000e+001,
                                        -2.95100000000000020e+001 );
    b[3] = lapack_make_complex_double( -1.47400000000000000e+001,
                                       -2.39999999999999990e+000 );
    b[11] = lapack_make_complex_double( 1.91700000000000020e+001,
                                        2.13299999999999980e+001 );
}

/* Auxiliary function: C interface to ztbtrs results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ztbtrs( lapack_complex_double *b, lapack_complex_double *b_i,
                           lapack_int info, lapack_int info_i, lapack_int ldb,
                           lapack_int nrhs )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldb*nrhs; i++ ) {
        failed += compare_complex_doubles(b[i],b_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
