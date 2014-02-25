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
* zpbtrf_1 is the test program for the C interface to LAPACK
* routine zpbtrf
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

static void init_scalars_zpbtrf( char *uplo, lapack_int *n, lapack_int *kd,
                                 lapack_int *ldab );
static void init_ab( lapack_int size, lapack_complex_double *ab );
static int compare_zpbtrf( lapack_complex_double *ab,
                           lapack_complex_double *ab_i, lapack_int info,
                           lapack_int info_i, lapack_int ldab, lapack_int n );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
    lapack_int n, n_i;
    lapack_int kd, kd_i;
    lapack_int ldab, ldab_i;
    lapack_int ldab_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_double *ab = NULL, *ab_i = NULL;
    lapack_complex_double *ab_save = NULL;
    lapack_complex_double *ab_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_zpbtrf( &uplo, &n, &kd, &ldab );
    ldab_r = n+2;
    uplo_i = uplo;
    n_i = n;
    kd_i = kd;
    ldab_i = ldab;

    /* Allocate memory for the LAPACK routine arrays */
    ab = (lapack_complex_double *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (lapack_complex_double *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the backup arrays */
    ab_save = (lapack_complex_double *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the row-major arrays */
    ab_r = (lapack_complex_double *)
        LAPACKE_malloc( (kd+1)*(n+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldab*n; i++ ) {
        ab_save[i] = ab[i];
    }

    /* Call the LAPACK routine */
    zpbtrf_( &uplo, &n, &kd, ab, &ldab, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    info_i = LAPACKE_zpbtrf_work( LAPACK_COL_MAJOR, uplo_i, n_i, kd_i, ab_i,
                                  ldab_i );

    failed = compare_zpbtrf( ab, ab_i, info, info_i, ldab, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zpbtrf\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zpbtrf\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    info_i = LAPACKE_zpbtrf( LAPACK_COL_MAJOR, uplo_i, n_i, kd_i, ab_i,
                             ldab_i );

    failed = compare_zpbtrf( ab, ab_i, info, info_i, ldab, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zpbtrf\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zpbtrf\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, kd+1, n, ab_i, ldab, ab_r, n+2 );
    info_i = LAPACKE_zpbtrf_work( LAPACK_ROW_MAJOR, uplo_i, n_i, kd_i, ab_r,
                                  ldab_r );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, kd+1, n, ab_r, n+2, ab_i, ldab );

    failed = compare_zpbtrf( ab, ab_i, info, info_i, ldab, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zpbtrf\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zpbtrf\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, kd+1, n, ab_i, ldab, ab_r, n+2 );
    info_i = LAPACKE_zpbtrf( LAPACK_ROW_MAJOR, uplo_i, n_i, kd_i, ab_r,
                             ldab_r );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, kd+1, n, ab_r, n+2, ab_i, ldab );

    failed = compare_zpbtrf( ab, ab_i, info, info_i, ldab, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zpbtrf\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zpbtrf\n" );
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
    if( ab_save != NULL ) {
        LAPACKE_free( ab_save );
    }

    return 0;
}

/* Auxiliary function: zpbtrf scalar parameters initialization */
static void init_scalars_zpbtrf( char *uplo, lapack_int *n, lapack_int *kd,
                                 lapack_int *ldab )
{
    *uplo = 'L';
    *n = 4;
    *kd = 1;
    *ldab = 9;

    return;
}

/* Auxiliary functions: zpbtrf array parameters initialization */
static void init_ab( lapack_int size, lapack_complex_double *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    ab[0] = lapack_make_complex_double( 9.39000000000000060e+000,
                                        0.00000000000000000e+000 );
    ab[9] = lapack_make_complex_double( 1.68999999999999990e+000,
                                        0.00000000000000000e+000 );
    ab[18] = lapack_make_complex_double( 2.64999999999999990e+000,
                                         0.00000000000000000e+000 );
    ab[27] = lapack_make_complex_double( 2.16999999999999990e+000,
                                         0.00000000000000000e+000 );
    ab[1] = lapack_make_complex_double( 1.08000000000000010e+000,
                                        1.73000000000000000e+000 );
    ab[10] = lapack_make_complex_double( -4.00000000000000010e-002,
                                         -2.89999999999999980e-001 );
    ab[19] = lapack_make_complex_double( -3.30000000000000020e-001,
                                         -2.24000000000000020e+000 );
    ab[28] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
}

/* Auxiliary function: C interface to zpbtrf results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zpbtrf( lapack_complex_double *ab,
                           lapack_complex_double *ab_i, lapack_int info,
                           lapack_int info_i, lapack_int ldab, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldab*n; i++ ) {
        failed += compare_complex_doubles(ab[i],ab_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
