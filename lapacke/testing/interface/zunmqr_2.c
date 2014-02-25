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
* zunmqr_2 is the test program for the C interface to LAPACK
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
    *m = 6;
    *n = 2;
    *k = 4;
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
    a[0] = lapack_make_complex_double( -3.08700502105195800e+000,
                                       0.00000000000000000e+000 );
    a[8] = lapack_make_complex_double( -4.88499367417976620e-001,
                                       -1.14168910512461390e+000 );
    a[16] = lapack_make_complex_double( 3.77356043173210700e-001,
                                        -1.24372975548049110e+000 );
    a[24] = lapack_make_complex_double( -8.55165437696762120e-001,
                                        -7.07319873181135650e-001 );
    a[1] = lapack_make_complex_double( -3.26978431123342020e-001,
                                       4.23806608064021150e-001 );
    a[9] = lapack_make_complex_double( 1.51631604729093160e+000,
                                       0.00000000000000000e+000 );
    a[17] = lapack_make_complex_double( 1.37305509662697790e+000,
                                        -8.17629335421158570e-001 );
    a[25] = lapack_make_complex_double( -2.50862720262847670e-001,
                                        8.20348604045144430e-001 );
    a[2] = lapack_make_complex_double( 1.69172476430465120e-001,
                                       -7.98047673307239890e-002 );
    a[10] = lapack_make_complex_double( -4.53710486149829670e-001,
                                        -6.49149959135295050e-003 );
    a[18] = lapack_make_complex_double( -2.17134536255717010e+000,
                                        0.00000000000000000e+000 );
    a[26] = lapack_make_complex_double( -2.27267620328359900e-001,
                                        -2.95731405907064150e-001 );
    a[3] = lapack_make_complex_double( -1.05973629513027950e-001,
                                       7.26861860966964010e-002 );
    a[11] = lapack_make_complex_double( -2.73407174139646880e-001,
                                        9.78078838704348470e-002 );
    a[19] = lapack_make_complex_double( -2.91822737804995960e-001,
                                        4.88808144155306160e-001 );
    a[27] = lapack_make_complex_double( -2.35337610655542080e+000,
                                        0.00000000000000000e+000 );
    a[4] = lapack_make_complex_double( 1.72939632545932170e-001,
                                       1.60632640429298590e-001 );
    a[12] = lapack_make_complex_double( -3.23630471463261900e-001,
                                        1.23000700219973860e-001 );
    a[20] = lapack_make_complex_double( 2.72768506164479230e-001,
                                        4.69769330690376700e-002 );
    a[28] = lapack_make_complex_double( 7.05422688603155710e-001,
                                        2.51508056610988850e-001 );
    a[5] = lapack_make_complex_double( 2.69899674468747240e-001,
                                       -1.51670836485297000e-002 );
    a[13] = lapack_make_complex_double( -1.64593543935458470e-001,
                                        3.38900720348261240e-001 );
    a[21] = lapack_make_complex_double( 5.34839525361778920e-001,
                                        3.98829067784022160e-001 );
    a[29] = lapack_make_complex_double( 2.70306990523076050e-001,
                                        -7.26878326406570770e-002 );
}
static void init_tau( lapack_int size, lapack_complex_double *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    tau[0] = lapack_make_complex_double( 1.31098102965600650e+000,
                                         -2.62390243772255500e-001 );
    tau[1] = lapack_make_complex_double( 1.10510398911057410e+000,
                                         -4.50362538745018080e-001 );
    tau[2] = lapack_make_complex_double( 1.04025187161551910e+000,
                                         2.12175810726109660e-001 );
    tau[3] = lapack_make_complex_double( 1.18595901116610980e+000,
                                         2.01183600330743640e-001 );
}
static void init_c( lapack_int size, lapack_complex_double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    c[0] = lapack_make_complex_double( -1.54000000000000000e+000,
                                       7.60000000000000010e-001 );
    c[8] = lapack_make_complex_double( 3.16999999999999990e+000,
                                       -2.08999999999999990e+000 );
    c[1] = lapack_make_complex_double( 1.20000000000000000e-001,
                                       -1.91999999999999990e+000 );
    c[9] = lapack_make_complex_double( -6.53000000000000020e+000,
                                       4.17999999999999970e+000 );
    c[2] = lapack_make_complex_double( -9.08000000000000010e+000,
                                       -4.30999999999999960e+000 );
    c[10] = lapack_make_complex_double( 7.28000000000000020e+000,
                                        7.29999999999999980e-001 );
    c[3] = lapack_make_complex_double( 7.49000000000000020e+000,
                                       3.64999999999999990e+000 );
    c[11] = lapack_make_complex_double( 9.10000000000000030e-001,
                                        -3.97000000000000020e+000 );
    c[4] = lapack_make_complex_double( -5.62999999999999990e+000,
                                       -2.12000000000000010e+000 );
    c[12] = lapack_make_complex_double( -5.46000000000000000e+000,
                                        -1.63999999999999990e+000 );
    c[5] = lapack_make_complex_double( 2.37000000000000010e+000,
                                       8.02999999999999940e+000 );
    c[13] = lapack_make_complex_double( -2.83999999999999990e+000,
                                        -5.86000000000000030e+000 );
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
