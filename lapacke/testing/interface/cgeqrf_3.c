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
* cgeqrf_3 is the test program for the C interface to LAPACK
* routine cgeqrf
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

static void init_scalars_cgeqrf( lapack_int *m, lapack_int *n, lapack_int *lda,
                                 lapack_int *lwork );
static void init_a( lapack_int size, lapack_complex_float *a );
static void init_tau( lapack_int size, lapack_complex_float *tau );
static void init_work( lapack_int size, lapack_complex_float *work );
static int compare_cgeqrf( lapack_complex_float *a, lapack_complex_float *a_i,
                           lapack_complex_float *tau,
                           lapack_complex_float *tau_i, lapack_int info,
                           lapack_int info_i, lapack_int lda, lapack_int m,
                           lapack_int n );

int main(void)
{
    /* Local scalars */
    lapack_int m, m_i;
    lapack_int n, n_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    lapack_int lwork, lwork_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *a = NULL, *a_i = NULL;
    lapack_complex_float *tau = NULL, *tau_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    lapack_complex_float *a_save = NULL;
    lapack_complex_float *tau_save = NULL;
    lapack_complex_float *a_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_cgeqrf( &m, &n, &lda, &lwork );
    lda_r = n+2;
    m_i = m;
    n_i = n;
    lda_i = lda;
    lwork_i = lwork;

    /* Allocate memory for the LAPACK routine arrays */
    a = (lapack_complex_float *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_float) );
    tau = (lapack_complex_float *)
        LAPACKE_malloc( MIN(m,n) * sizeof(lapack_complex_float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_float *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_float) );
    tau_i = (lapack_complex_float *)
        LAPACKE_malloc( MIN(m,n) * sizeof(lapack_complex_float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    a_save = (lapack_complex_float *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_float) );
    tau_save = (lapack_complex_float *)
        LAPACKE_malloc( MIN(m,n) * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_float *)
        LAPACKE_malloc( m*(n+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_a( lda*n, a );
    init_tau( (MIN(m,n)), tau );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < lda*n; i++ ) {
        a_save[i] = a[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_save[i] = tau[i];
    }

    /* Call the LAPACK routine */
    cgeqrf_( &m, &n, a, &lda, tau, work, &lwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_i[i] = tau_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_cgeqrf_work( LAPACK_COL_MAJOR, m_i, n_i, a_i, lda_i, tau_i,
                                  work_i, lwork_i );

    failed = compare_cgeqrf( a, a_i, tau, tau_i, info, info_i, lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to cgeqrf\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to cgeqrf\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_i[i] = tau_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_cgeqrf( LAPACK_COL_MAJOR, m_i, n_i, a_i, lda_i, tau_i );

    failed = compare_cgeqrf( a, a_i, tau, tau_i, info, info_i, lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to cgeqrf\n" );
    } else {
        printf( "FAILED: column-major high-level interface to cgeqrf\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_i[i] = tau_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_cgeqrf_work( LAPACK_ROW_MAJOR, m_i, n_i, a_r, lda_r, tau_i,
                                  work_i, lwork_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, m, n, a_r, n+2, a_i, lda );

    failed = compare_cgeqrf( a, a_i, tau, tau_i, info, info_i, lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to cgeqrf\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to cgeqrf\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_i[i] = tau_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_cgeqrf( LAPACK_ROW_MAJOR, m_i, n_i, a_r, lda_r, tau_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, m, n, a_r, n+2, a_i, lda );

    failed = compare_cgeqrf( a, a_i, tau, tau_i, info, info_i, lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to cgeqrf\n" );
    } else {
        printf( "FAILED: row-major high-level interface to cgeqrf\n" );
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
    if( a_save != NULL ) {
        LAPACKE_free( a_save );
    }
    if( tau != NULL ) {
        LAPACKE_free( tau );
    }
    if( tau_i != NULL ) {
        LAPACKE_free( tau_i );
    }
    if( tau_save != NULL ) {
        LAPACKE_free( tau_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: cgeqrf scalar parameters initialization */
static void init_scalars_cgeqrf( lapack_int *m, lapack_int *n, lapack_int *lda,
                                 lapack_int *lwork )
{
    *m = 6;
    *n = 4;
    *lda = 8;
    *lwork = 1024;

    return;
}

/* Auxiliary functions: cgeqrf array parameters initialization */
static void init_a( lapack_int size, lapack_complex_float *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    a[0] = lapack_make_complex_float( 9.599999785e-001, -8.100000024e-001 );
    a[8] = lapack_make_complex_float( -2.999999933e-002, 9.599999785e-001 );
    a[16] = lapack_make_complex_float( -9.100000262e-001, 2.059999943e+000 );
    a[24] = lapack_make_complex_float( -5.000000075e-002, 4.099999964e-001 );
    a[1] = lapack_make_complex_float( -9.800000191e-001, 1.980000019e+000 );
    a[9] = lapack_make_complex_float( -1.200000048e+000, 1.899999976e-001 );
    a[17] = lapack_make_complex_float( -6.600000262e-001, 4.199999869e-001 );
    a[25] = lapack_make_complex_float( -8.100000024e-001, 5.600000024e-001 );
    a[2] = lapack_make_complex_float( 6.200000048e-001, -4.600000083e-001 );
    a[10] = lapack_make_complex_float( 1.009999990e+000, 1.999999955e-002 );
    a[18] = lapack_make_complex_float( 6.299999952e-001, -1.700000018e-001 );
    a[26] = lapack_make_complex_float( -1.110000014e+000, 6.000000238e-001 );
    a[3] = lapack_make_complex_float( -3.700000048e-001, 3.799999952e-001 );
    a[11] = lapack_make_complex_float( 1.899999976e-001, -5.400000215e-001 );
    a[19] = lapack_make_complex_float( -9.800000191e-001, -3.600000143e-001 );
    a[27] = lapack_make_complex_float( 2.199999988e-001, -2.000000030e-001 );
    a[4] = lapack_make_complex_float( 8.299999833e-001, 5.099999905e-001 );
    a[12] = lapack_make_complex_float( 2.000000030e-001, 9.999999776e-003 );
    a[20] = lapack_make_complex_float( -1.700000018e-001, -4.600000083e-001 );
    a[28] = lapack_make_complex_float( 1.470000029e+000, 1.590000033e+000 );
    a[5] = lapack_make_complex_float( 1.080000043e+000, -2.800000012e-001 );
    a[13] = lapack_make_complex_float( 2.000000030e-001, -1.199999973e-001 );
    a[21] = lapack_make_complex_float( -7.000000030e-002, 1.230000019e+000 );
    a[29] = lapack_make_complex_float( 2.599999905e-001, 2.599999905e-001 );
}
static void init_tau( lapack_int size, lapack_complex_float *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}

/* Auxiliary function: C interface to cgeqrf results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_cgeqrf( lapack_complex_float *a, lapack_complex_float *a_i,
                           lapack_complex_float *tau,
                           lapack_complex_float *tau_i, lapack_int info,
                           lapack_int info_i, lapack_int lda, lapack_int m,
                           lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < lda*n; i++ ) {
        failed += compare_complex_floats(a[i],a_i[i]);
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        failed += compare_complex_floats(tau[i],tau_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
