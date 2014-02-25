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
* cunmlq_1 is the test program for the C interface to LAPACK
* routine cunmlq
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

static void init_scalars_cunmlq( char *side, char *trans, lapack_int *m,
                                 lapack_int *n, lapack_int *k, lapack_int *lda,
                                 lapack_int *ldc, lapack_int *lwork );
static void init_a( lapack_int size, lapack_complex_float *a );
static void init_tau( lapack_int size, lapack_complex_float *tau );
static void init_c( lapack_int size, lapack_complex_float *c );
static void init_work( lapack_int size, lapack_complex_float *work );
static int compare_cunmlq( lapack_complex_float *c, lapack_complex_float *c_i,
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
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *a = NULL, *a_i = NULL;
    lapack_complex_float *tau = NULL, *tau_i = NULL;
    lapack_complex_float *c = NULL, *c_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    lapack_complex_float *c_save = NULL;
    lapack_complex_float *a_r = NULL;
    lapack_complex_float *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_cunmlq( &side, &trans, &m, &n, &k, &lda, &ldc, &lwork );
    lda_r = m+2;
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
    a = (lapack_complex_float *)
        LAPACKE_malloc( lda*m * sizeof(lapack_complex_float) );
    tau = (lapack_complex_float *)
        LAPACKE_malloc( k * sizeof(lapack_complex_float) );
    c = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_float *)
        LAPACKE_malloc( lda*m * sizeof(lapack_complex_float) );
    tau_i = (lapack_complex_float *)
        LAPACKE_malloc( k * sizeof(lapack_complex_float) );
    c_i = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    c_save = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_float *)
        LAPACKE_malloc( k*(m+2) * sizeof(lapack_complex_float) );
    c_r = (lapack_complex_float *)
        LAPACKE_malloc( m*(n+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_a( lda*m, a );
    init_tau( k, tau );
    init_c( ldc*n, c );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldc*n; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    cunmlq_( &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork,
             &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*m; i++ ) {
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
    info_i = LAPACKE_cunmlq_work( LAPACK_COL_MAJOR, side_i, trans_i, m_i, n_i,
                                  k_i, a_i, lda_i, tau_i, c_i, ldc_i, work_i,
                                  lwork_i );

    failed = compare_cunmlq( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to cunmlq\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to cunmlq\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*m; i++ ) {
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
    info_i = LAPACKE_cunmlq( LAPACK_COL_MAJOR, side_i, trans_i, m_i, n_i, k_i,
                             a_i, lda_i, tau_i, c_i, ldc_i );

    failed = compare_cunmlq( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to cunmlq\n" );
    } else {
        printf( "FAILED: column-major high-level interface to cunmlq\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*m; i++ ) {
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

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, k, m, a_i, lda, a_r, m+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_cunmlq_work( LAPACK_ROW_MAJOR, side_i, trans_i, m_i, n_i,
                                  k_i, a_r, lda_r, tau_i, c_r, ldc_r, work_i,
                                  lwork_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_cunmlq( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to cunmlq\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to cunmlq\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*m; i++ ) {
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
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, k, m, a_i, lda, a_r, m+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_cunmlq( LAPACK_ROW_MAJOR, side_i, trans_i, m_i, n_i, k_i,
                             a_r, lda_r, tau_i, c_r, ldc_r );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_cunmlq( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to cunmlq\n" );
    } else {
        printf( "FAILED: row-major high-level interface to cunmlq\n" );
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

/* Auxiliary function: cunmlq scalar parameters initialization */
static void init_scalars_cunmlq( char *side, char *trans, lapack_int *m,
                                 lapack_int *n, lapack_int *k, lapack_int *lda,
                                 lapack_int *ldc, lapack_int *lwork )
{
    *side = 'L';
    *trans = 'C';
    *m = 4;
    *n = 2;
    *k = 3;
    *lda = 8;
    *ldc = 8;
    *lwork = 512;

    return;
}

/* Auxiliary functions: cunmlq array parameters initialization */
static void init_a( lapack_int size, lapack_complex_float *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    a[0] = lapack_make_complex_float( -2.225511312e+000, 0.000000000e+000 );
    a[8] = lapack_make_complex_float( 2.438442558e-001, -3.082070053e-001 );
    a[16] = lapack_make_complex_float( -2.741365135e-001, -2.309664786e-001 );
    a[24] = lapack_make_complex_float( 5.807709098e-001, 3.468663692e-001 );
    a[1] = lapack_make_complex_float( 8.207554817e-001, 1.238457084e+000 );
    a[9] = lapack_make_complex_float( 1.688100934e+000, 0.000000000e+000 );
    a[17] = lapack_make_complex_float( -1.936414540e-001, 5.429518223e-001 );
    a[25] = lapack_make_complex_float( 2.789084315e-001, -2.203175873e-001 );
    a[2] = lapack_make_complex_float( 1.033455133e-003, -6.822252274e-001 );
    a[10] = lapack_make_complex_float( 7.747511864e-001, -6.154725552e-001 );
    a[18] = lapack_make_complex_float( -1.590258121e+000, 0.000000000e+000 );
    a[26] = lapack_make_complex_float( -1.267668605e-001, 1.109846011e-001 );
}
static void init_tau( lapack_int size, lapack_complex_float *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    tau[0] = lapack_make_complex_float( 1.125813842e+000, 1.617605835e-001 );
    tau[1] = lapack_make_complex_float( 1.099053621e+000, 5.468589664e-001 );
    tau[2] = lapack_make_complex_float( 1.132925749e+000, -9.590539932e-001 );
}
static void init_c( lapack_int size, lapack_complex_float *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    c[0] = lapack_make_complex_float( 6.066021919e-001, -8.537364006e-002 );
    c[8] = lapack_make_complex_float( -2.170287848e+000, 1.199724317e+000 );
    c[1] = lapack_make_complex_float( 5.216747761e+000, -2.512397289e+000 );
    c[9] = lapack_make_complex_float( -2.377178669e+000, 2.987456560e+000 );
    c[2] = lapack_make_complex_float( 6.293161869e+000, -7.861097813e+000 );
    c[10] = lapack_make_complex_float( 1.214990690e-001, 4.587235153e-001 );
    c[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    c[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}

/* Auxiliary function: C interface to cunmlq results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_cunmlq( lapack_complex_float *c, lapack_complex_float *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldc*n; i++ ) {
        failed += compare_complex_floats(c[i],c_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
