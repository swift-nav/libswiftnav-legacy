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
* zgecon_1 is the test program for the C interface to LAPACK
* routine zgecon
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

static void init_scalars_zgecon( char *norm, lapack_int *n, lapack_int *lda,
                                 double *anorm );
static void init_a( lapack_int size, lapack_complex_double *a );
static void init_work( lapack_int size, lapack_complex_double *work );
static void init_rwork( lapack_int size, double *rwork );
static int compare_zgecon( double rcond, double rcond_i, lapack_int info,
                           lapack_int info_i );

int main(void)
{
    /* Local scalars */
    char norm, norm_i;
    lapack_int n, n_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    double anorm, anorm_i;
    double rcond, rcond_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_double *a = NULL, *a_i = NULL;
    lapack_complex_double *work = NULL, *work_i = NULL;
    double *rwork = NULL, *rwork_i = NULL;
    lapack_complex_double *a_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_zgecon( &norm, &n, &lda, &anorm );
    lda_r = n+2;
    norm_i = norm;
    n_i = n;
    lda_i = lda;
    anorm_i = anorm;

    /* Allocate memory for the LAPACK routine arrays */
    a = (lapack_complex_double *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_double) );
    work = (lapack_complex_double *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_double) );
    rwork = (double *)LAPACKE_malloc( 2*n * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_double *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_double) );
    work_i = (lapack_complex_double *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_double) );
    rwork_i = (double *)LAPACKE_malloc( 2*n * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_a( lda*n, a );
    init_work( 2*n, work );
    init_rwork( 2*n, rwork );

    /* Call the LAPACK routine */
    zgecon_( &norm, &n, a, &lda, &anorm, &rcond, work, rwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_zgecon_work( LAPACK_COL_MAJOR, norm_i, n_i, a_i, lda_i,
                                  anorm_i, &rcond_i, work_i, rwork_i );

    failed = compare_zgecon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zgecon\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zgecon\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_zgecon( LAPACK_COL_MAJOR, norm_i, n_i, a_i, lda_i, anorm_i,
                             &rcond_i );

    failed = compare_zgecon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zgecon\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zgecon\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        rwork_i[i] = rwork[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_zgecon_work( LAPACK_ROW_MAJOR, norm_i, n_i, a_r, lda_r,
                                  anorm_i, &rcond_i, work_i, rwork_i );

    failed = compare_zgecon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zgecon\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zgecon\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        rwork_i[i] = rwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_zgecon( LAPACK_ROW_MAJOR, norm_i, n_i, a_r, lda_r, anorm_i,
                             &rcond_i );

    failed = compare_zgecon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zgecon\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zgecon\n" );
    }

    /* Release memory */
    if( a != NULL ) {
        LAPACKE_free( a );
    }
    if( a_i != NULL ) {
        LAPACKE_free( a_i );
    }
    if( a_r != NULL ) {
        LAPACKE_free( a_r );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }
    if( rwork != NULL ) {
        LAPACKE_free( rwork );
    }
    if( rwork_i != NULL ) {
        LAPACKE_free( rwork_i );
    }

    return 0;
}

/* Auxiliary function: zgecon scalar parameters initialization */
static void init_scalars_zgecon( char *norm, lapack_int *n, lapack_int *lda,
                                 double *anorm )
{
    *norm = '1';
    *n = 4;
    *lda = 8;
    *anorm = 1.28838218636040820e+001;

    return;
}

/* Auxiliary functions: zgecon array parameters initialization */
static void init_a( lapack_int size, lapack_complex_double *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    a[0] = lapack_make_complex_double( -3.29000000000000000e+000,
                                       -2.39000000000000010e+000 );
    a[8] = lapack_make_complex_double( -1.90999999999999990e+000,
                                       4.41999999999999990e+000 );
    a[16] = lapack_make_complex_double( -1.40000000000000010e-001,
                                        -1.35000000000000010e+000 );
    a[24] = lapack_make_complex_double( 1.72000000000000000e+000,
                                        1.35000000000000010e+000 );
    a[1] = lapack_make_complex_double( 2.37612026946940610e-001,
                                       2.55959652157085600e-001 );
    a[9] = lapack_make_complex_double( 4.89518063400297530e+000,
                                       -7.11362223485444090e-001 );
    a[17] = lapack_make_complex_double( -4.62279846639493950e-001,
                                        1.69661058768036190e+000 );
    a[25] = lapack_make_complex_double( 1.22685284406332770e+000,
                                        6.18973161911442920e-001 );
    a[2] = lapack_make_complex_double( -1.01952080889200600e-001,
                                       -7.01013533943711350e-001 );
    a[10] = lapack_make_complex_double( -6.69149643194321130e-001,
                                        3.68869854797165500e-001 );
    a[18] = lapack_make_complex_double( -5.14141091381023150e+000,
                                        -1.12996973466095300e+000 );
    a[26] = lapack_make_complex_double( 9.98257991519944880e-001,
                                        3.85015227576377740e-001 );
    a[3] = lapack_make_complex_double( -5.35854670359574790e-001,
                                       2.70727252935982880e-001 );
    a[11] = lapack_make_complex_double( -2.04021779458707920e-001,
                                        8.60118067998404510e-001 );
    a[19] = lapack_make_complex_double( 8.23304706394011040e-003,
                                        1.21063681991063470e-001 );
    a[27] = lapack_make_complex_double( 1.48239181074841820e-001,
                                        -1.25223998607584600e-001 );
}
static void init_work( lapack_int size, lapack_complex_double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
}
static void init_rwork( lapack_int size, double *rwork ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        rwork[i] = 0;
    }
}

/* Auxiliary function: C interface to zgecon results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zgecon( double rcond, double rcond_i, lapack_int info,
                           lapack_int info_i )
{
    int failed = 0;
    failed += compare_doubles(rcond,rcond_i);
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
