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
* zunmqr_1 is the test program for the C interface to LAPACK
* routine zunmqr
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

static void init_scalars_zunmqr( char *side, char *trans, lapack_int *m,
                                 lapack_int *n, lapack_int *k, lapack_int *lda,
                                 lapack_int *ldc, lapack_int *lwork );
static void init_a( lapack_int size, lapack_complex_double *a );
static void init_tau( lapack_int size, lapack_complex_double *tau );
static void init_c( lapack_int size, lapack_complex_double *c );
static void init_work( lapack_int size, lapack_complex_double *work );
static int compare_zunmqr( lapack_complex_double *c, lapack_complex_double *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int n );

int main(void)
{
    /* Local scalars */
    char side, side_i;
    char trans, trans_i;
    lapack_int m, m_i;
    lapack_int n, n_i;
    lapack_int k, k_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    lapack_int ldc, ldc_i;
    lapack_int ldc_r;
    lapack_int lwork, lwork_i;
    lapack_int info, info_i;
    /* Declare scalars */
    lapack_int r;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_double *a = NULL, *a_i = NULL;
    lapack_complex_double *tau = NULL, *tau_i = NULL;
    lapack_complex_double *c = NULL, *c_i = NULL;
    lapack_complex_double *work = NULL, *work_i = NULL;
    lapack_complex_double *c_save = NULL;
    lapack_complex_double *a_r = NULL;
    lapack_complex_double *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_zunmqr( &side, &trans, &m, &n, &k, &lda, &ldc, &lwork );
    r = LAPACKE_lsame( side, 'l' ) ? m : n;
    lda_r = k+2;
    ldc_r = n+2;
    side_i = side;
    trans_i = trans;
    m_i = m;
    n_i = n;
    k_i = k;
    lda_i = lda;
    ldc_i = ldc;
    lwork_i = lwork;

    /* Allocate memory for the LAPACK routine arrays */
    a = (lapack_complex_double *)
        LAPACKE_malloc( lda*k * sizeof(lapack_complex_double) );
    tau = (lapack_complex_double *)
        LAPACKE_malloc( k * sizeof(lapack_complex_double) );
    c = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );
    work = (lapack_complex_double *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_double) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_double *)
        LAPACKE_malloc( lda*k * sizeof(lapack_complex_double) );
    tau_i = (lapack_complex_double *)
        LAPACKE_malloc( k * sizeof(lapack_complex_double) );
    c_i = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );
    work_i = (lapack_complex_double *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_double) );

    /* Allocate memory for the backup arrays */
    c_save = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_double *)
        LAPACKE_malloc( r*(k+2) * sizeof(lapack_complex_double) );
    c_r = (lapack_complex_double *)
        LAPACKE_malloc( m*(n+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_a( lda*k, a );
    init_tau( k, tau );
    init_c( ldc*n, c );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldc*n; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    zunmqr_( &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork,
             &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*k; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < k; i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_zunmqr_work( LAPACK_COL_MAJOR, side_i, trans_i, m_i, n_i,
                                  k_i, a_i, lda_i, tau_i, c_i, ldc_i, work_i,
                                  lwork_i );

    failed = compare_zunmqr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zunmqr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zunmqr\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*k; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < k; i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_zunmqr( LAPACK_COL_MAJOR, side_i, trans_i, m_i, n_i, k_i,
                             a_i, lda_i, tau_i, c_i, ldc_i );

    failed = compare_zunmqr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zunmqr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zunmqr\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*k; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < k; i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, r, k, a_i, lda, a_r, k+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_zunmqr_work( LAPACK_ROW_MAJOR, side_i, trans_i, m_i, n_i,
                                  k_i, a_r, lda_r, tau_i, c_r, ldc_r, work_i,
                                  lwork_i );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_zunmqr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zunmqr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zunmqr\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*k; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < k; i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, r, k, a_i, lda, a_r, k+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_zunmqr( LAPACK_ROW_MAJOR, side_i, trans_i, m_i, n_i, k_i,
                             a_r, lda_r, tau_i, c_r, ldc_r );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_zunmqr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zunmqr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zunmqr\n" );
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
    if( tau != NULL ) {
        LAPACKE_free( tau );
    }
    if( tau_i != NULL ) {
        LAPACKE_free( tau_i );
    }
    if( c != NULL ) {
        LAPACKE_free( c );
    }
    if( c_i != NULL ) {
        LAPACKE_free( c_i );
    }
    if( c_r != NULL ) {
        LAPACKE_free( c_r );
    }
    if( c_save != NULL ) {
        LAPACKE_free( c_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: zunmqr scalar parameters initialization */
static void init_scalars_zunmqr( char *side, char *trans, lapack_int *m,
                                 lapack_int *n, lapack_int *k, lapack_int *lda,
                                 lapack_int *ldc, lapack_int *lwork )
{
    *side = 'L';
    *trans = 'C';
    *m = 5;
    *n = 2;
    *k = 3;
    *lda = 8;
    *ldc = 8;
    *lwork = 512;

    return;
}

/* Auxiliary functions: zunmqr array parameters initialization */
static void init_a( lapack_int size, lapack_complex_double *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    a[0] = lapack_make_complex_double( -2.28700677742764880e+000,
                                       0.00000000000000000e+000 );
    a[8] = lapack_make_complex_double( -7.23871925671358720e-001,
                                       5.60514298712240460e-001 );
    a[16] = lapack_make_complex_double( 7.56971958757025500e-001,
                                        4.51944441180257450e-001 );
    a[1] = lapack_make_complex_double( -1.41987831087444230e-001,
                                       8.15431066587506830e-003 );
    a[9] = lapack_make_complex_double( 1.89500215202093210e+000,
                                       0.00000000000000000e+000 );
    a[17] = lapack_make_complex_double( -3.62042107080837280e-001,
                                        -4.78657240305253420e-001 );
    a[2] = lapack_make_complex_double( -9.06131808634396410e-002,
                                       -5.90946044610632520e-002 );
    a[10] = lapack_make_complex_double( -4.73439672577799340e-001,
                                        1.68439572892672610e-001 );
    a[18] = lapack_make_complex_double( 9.39282936471031560e-001,
                                        0.00000000000000000e+000 );
    a[3] = lapack_make_complex_double( 4.08210340556465660e-001,
                                       -1.17792243052802140e-001 );
    a[11] = lapack_make_complex_double( 2.58521862467017920e-001,
                                        2.66720382327443330e-002 );
    a[19] = lapack_make_complex_double( 4.50972094106415690e-002,
                                        5.81513806297987880e-001 );
    a[4] = lapack_make_complex_double( -8.95656153586470990e-002,
                                       3.39695746702570130e-001 );
    a[12] = lapack_make_complex_double( -6.87183019669032950e-001,
                                        3.17086911183261100e-002 );
    a[20] = lapack_make_complex_double( -1.17064906202052320e-001,
                                        -1.20657309265043160e-002 );
}
static void init_tau( lapack_int size, lapack_complex_double *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    tau[0] = lapack_make_complex_double( 1.34980219905592680e+000,
                                         -4.45997803796306660e-001 );
    tau[1] = lapack_make_complex_double( 1.04037745348862610e+000,
                                         2.79135088926636740e-001 );
    tau[2] = lapack_make_complex_double( 1.47028895497071630e+000,
                                         9.97710593712547170e-002 );
}
static void init_c( lapack_int size, lapack_complex_double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    c[0] = lapack_make_complex_double( -8.49999999999999980e-001,
                                       -1.62999999999999990e+000 );
    c[8] = lapack_make_complex_double( 2.49000000000000020e+000,
                                       4.00999999999999980e+000 );
    c[1] = lapack_make_complex_double( -2.16000000000000010e+000,
                                       3.52000000000000000e+000 );
    c[9] = lapack_make_complex_double( -1.40000000000000010e-001,
                                       7.98000000000000040e+000 );
    c[2] = lapack_make_complex_double( 4.57000000000000030e+000,
                                       -5.71000000000000000e+000 );
    c[10] = lapack_make_complex_double( 8.35999999999999940e+000,
                                        -2.80000000000000030e-001 );
    c[3] = lapack_make_complex_double( 6.37999999999999990e+000,
                                       -7.40000000000000040e+000 );
    c[11] = lapack_make_complex_double( -3.54999999999999980e+000,
                                        1.29000000000000000e+000 );
    c[4] = lapack_make_complex_double( 8.41000000000000010e+000,
                                       9.39000000000000060e+000 );
    c[12] = lapack_make_complex_double( -6.71999999999999980e+000,
                                        5.03000000000000020e+000 );
}
static void init_work( lapack_int size, lapack_complex_double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
}

/* Auxiliary function: C interface to zunmqr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zunmqr( lapack_complex_double *c, lapack_complex_double *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldc*n; i++ ) {
        failed += compare_complex_doubles(c[i],c_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
