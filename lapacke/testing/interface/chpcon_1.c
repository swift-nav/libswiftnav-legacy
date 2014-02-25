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
* chpcon_1 is the test program for the C interface to LAPACK
* routine chpcon
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

static void init_scalars_chpcon( char *uplo, lapack_int *n, float *anorm );
static void init_ap( lapack_int size, lapack_complex_float *ap );
static void init_ipiv( lapack_int size, lapack_int *ipiv );
static void init_work( lapack_int size, lapack_complex_float *work );
static int compare_chpcon( float rcond, float rcond_i, lapack_int info,
                           lapack_int info_i );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
    lapack_int n, n_i;
    float anorm, anorm_i;
    float rcond, rcond_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *ap = NULL, *ap_i = NULL;
    lapack_int *ipiv = NULL, *ipiv_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    lapack_complex_float *ap_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_chpcon( &uplo, &n, &anorm );
    uplo_i = uplo;
    n_i = n;
    anorm_i = anorm;

    /* Allocate memory for the LAPACK routine arrays */
    ap = (lapack_complex_float *)
        LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(lapack_complex_float) );
    ipiv = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    ap_i = (lapack_complex_float *)
        LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(lapack_complex_float) );
    ipiv_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    ap_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+1)/2 * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_ap( (n*(n+1)/2), ap );
    init_ipiv( n, ipiv );
    init_work( 2*n, work );

    /* Call the LAPACK routine */
    chpcon_( &uplo, &n, ap, ipiv, &anorm, &rcond, work, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_chpcon_work( LAPACK_COL_MAJOR, uplo_i, n_i, ap_i, ipiv_i,
                                  anorm_i, &rcond_i, work_i );

    failed = compare_chpcon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to chpcon\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to chpcon\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_chpcon( LAPACK_COL_MAJOR, uplo_i, n_i, ap_i, ipiv_i,
                             anorm_i, &rcond_i );

    failed = compare_chpcon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to chpcon\n" );
    } else {
        printf( "FAILED: column-major high-level interface to chpcon\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_cpp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    info_i = LAPACKE_chpcon_work( LAPACK_ROW_MAJOR, uplo_i, n_i, ap_r, ipiv_i,
                                  anorm_i, &rcond_i, work_i );

    failed = compare_chpcon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to chpcon\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to chpcon\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_cpp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    info_i = LAPACKE_chpcon( LAPACK_ROW_MAJOR, uplo_i, n_i, ap_r, ipiv_i,
                             anorm_i, &rcond_i );

    failed = compare_chpcon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to chpcon\n" );
    } else {
        printf( "FAILED: row-major high-level interface to chpcon\n" );
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
    if( ipiv != NULL ) {
        LAPACKE_free( ipiv );
    }
    if( ipiv_i != NULL ) {
        LAPACKE_free( ipiv_i );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: chpcon scalar parameters initialization */
static void init_scalars_chpcon( char *uplo, lapack_int *n, float *anorm )
{
    *uplo = 'L';
    *n = 4;
    *anorm = 1.466419792e+001;

    return;
}

/* Auxiliary functions: chpcon array parameters initialization */
static void init_ap( lapack_int size, lapack_complex_float *ap ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ap[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    ap[0] = lapack_make_complex_float( -1.360000014e+000, 0.000000000e+000 );
    ap[1] = lapack_make_complex_float( 3.910000086e+000, -1.500000000e+000 );
    ap[2] = lapack_make_complex_float( 3.100287914e-001, 4.333020747e-002 );
    ap[3] = lapack_make_complex_float( -1.518120319e-001, 3.742958307e-001 );
    ap[4] = lapack_make_complex_float( -1.840000033e+000, 0.000000000e+000 );
    ap[5] = lapack_make_complex_float( 5.637050867e-001, 2.850349545e-001 );
    ap[6] = lapack_make_complex_float( 3.396582901e-001, 3.031452000e-002 );
    ap[7] = lapack_make_complex_float( -5.417624474e+000, 0.000000000e+000 );
    ap[8] = lapack_make_complex_float( 2.997244298e-001, 1.578268558e-001 );
    ap[9] = lapack_make_complex_float( -7.102810383e+000, 0.000000000e+000 );
}
static void init_ipiv( lapack_int size, lapack_int *ipiv ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ipiv[i] = 0;
    }
    ipiv[0] = -4;
    ipiv[1] = -4;
    ipiv[2] = 3;
    ipiv[3] = 4;
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}

/* Auxiliary function: C interface to chpcon results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_chpcon( float rcond, float rcond_i, lapack_int info,
                           lapack_int info_i )
{
    int failed = 0;
    failed += compare_floats(rcond,rcond_i);
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
