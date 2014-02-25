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
* zunmbr_1 is the test program for the C interface to LAPACK
* routine zunmbr
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

static void init_scalars_zunmbr( char *vect, char *side, char *trans,
                                 lapack_int *m, lapack_int *n, lapack_int *k,
                                 lapack_int *lda, lapack_int *ldc,
                                 lapack_int *lwork );
static void init_a( lapack_int size, lapack_complex_double *a );
static void init_tau( lapack_int size, lapack_complex_double *tau );
static void init_c( lapack_int size, lapack_complex_double *c );
static void init_work( lapack_int size, lapack_complex_double *work );
static int compare_zunmbr( lapack_complex_double *c, lapack_complex_double *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int n );

int main(void)
{
    /* Local scalars */
    char vect, vect_i;
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
    lapack_int nq;
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
    init_scalars_zunmbr( &vect, &side, &trans, &m, &n, &k, &lda, &ldc, &lwork );
    nq = LAPACKE_lsame( side, 'l' ) ? m : n;
    r = LAPACKE_lsame( vect, 'q' ) ? nq : MIN(nq,k);
    lda_r = MIN(nq,k)+2;
    ldc_r = n+2;
    vect_i = vect;
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
        LAPACKE_malloc( (lda*(MIN(nq,k))) * sizeof(lapack_complex_double) );
    tau = (lapack_complex_double *)
        LAPACKE_malloc( MIN(nq,k) * sizeof(lapack_complex_double) );
    c = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );
    work = (lapack_complex_double *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_double) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_double *)
        LAPACKE_malloc( (lda*(MIN(nq,k))) * sizeof(lapack_complex_double) );
    tau_i = (lapack_complex_double *)
        LAPACKE_malloc( MIN(nq,k) * sizeof(lapack_complex_double) );
    c_i = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );
    work_i = (lapack_complex_double *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_double) );

    /* Allocate memory for the backup arrays */
    c_save = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_double *)
        LAPACKE_malloc( (r*(MIN(nq,k)+2)) * sizeof(lapack_complex_double) );
    c_r = (lapack_complex_double *)
        LAPACKE_malloc( m*(n+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_a( lda*(MIN(nq,k)), a );
    init_tau( (MIN(nq,k)), tau );
    init_c( ldc*n, c );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldc*n; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    zunmbr_( &vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work,
             &lwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*(MIN(nq,k)); i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < (MIN(nq,k)); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_zunmbr_work( LAPACK_COL_MAJOR, vect_i, side_i, trans_i,
                                  m_i, n_i, k_i, a_i, lda_i, tau_i, c_i, ldc_i,
                                  work_i, lwork_i );

    failed = compare_zunmbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zunmbr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zunmbr\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*(MIN(nq,k)); i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < (MIN(nq,k)); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_zunmbr( LAPACK_COL_MAJOR, vect_i, side_i, trans_i, m_i,
                             n_i, k_i, a_i, lda_i, tau_i, c_i, ldc_i );

    failed = compare_zunmbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zunmbr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zunmbr\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*(MIN(nq,k)); i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < (MIN(nq,k)); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, r, MIN(nq, k ), a_i, lda, a_r, MIN(nq,
                       k)+2);
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_zunmbr_work( LAPACK_ROW_MAJOR, vect_i, side_i, trans_i,
                                  m_i, n_i, k_i, a_r, lda_r, tau_i, c_r, ldc_r,
                                  work_i, lwork_i );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_zunmbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zunmbr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zunmbr\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*(MIN(nq,k)); i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < (MIN(nq,k)); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, r, MIN(nq, k ), a_i, lda, a_r, MIN(nq,
                       k)+2);
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_zunmbr( LAPACK_ROW_MAJOR, vect_i, side_i, trans_i, m_i,
                             n_i, k_i, a_r, lda_r, tau_i, c_r, ldc_r );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_zunmbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zunmbr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zunmbr\n" );
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

/* Auxiliary function: zunmbr scalar parameters initialization */
static void init_scalars_zunmbr( char *vect, char *side, char *trans,
                                 lapack_int *m, lapack_int *n, lapack_int *k,
                                 lapack_int *lda, lapack_int *ldc,
                                 lapack_int *lwork )
{
    *vect = 'Q';
    *side = 'R';
    *trans = 'N';
    *m = 6;
    *n = 4;
    *k = 4;
    *lda = 8;
    *ldc = 8;
    *lwork = 1024;

    return;
}

/* Auxiliary functions: zunmbr array parameters initialization */
static void init_a( lapack_int size, lapack_complex_double *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    a[0] = lapack_make_complex_double( -3.08700502105195800e+000,
                                       0.00000000000000000e+000 );
    a[8] = lapack_make_complex_double( 2.11257100745583950e+000,
                                       0.00000000000000000e+000 );
    a[16] = lapack_make_complex_double( 5.43341107944031730e-002,
                                        4.54311849677352220e-001 );
    a[24] = lapack_make_complex_double( 3.75743827925403060e-001,
                                        1.07008730409452330e-001 );
    a[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    a[9] = lapack_make_complex_double( -2.06603927667906810e+000,
                                       0.00000000000000000e+000 );
    a[17] = lapack_make_complex_double( -1.26281010665522380e+000,
                                        0.00000000000000000e+000 );
    a[25] = lapack_make_complex_double( 2.82771782873273430e-002,
                                        1.65005610304937070e-001 );
    a[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    a[10] = lapack_make_complex_double( -2.80478799113691670e-001,
                                        -4.12446107471391210e-001 );
    a[18] = lapack_make_complex_double( -1.87312889112571160e+000,
                                        0.00000000000000000e+000 );
    a[26] = lapack_make_complex_double( 1.61263387280039240e+000,
                                        0.00000000000000000e+000 );
    a[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    a[11] = lapack_make_complex_double( 2.10347237263873110e-001,
                                        -4.46076099427661470e-001 );
    a[19] = lapack_make_complex_double( -5.70841942484137550e-001,
                                        6.43744629522165630e-002 );
    a[27] = lapack_make_complex_double( -2.00218286620699090e+000,
                                        0.00000000000000000e+000 );
}
static void init_tau( lapack_int size, lapack_complex_double *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    tau[0] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    tau[1] = lapack_make_complex_double( 1.09819823812611190e+000,
                                         5.15816216039656440e-001 );
    tau[2] = lapack_make_complex_double( 1.45515808833705250e+000,
                                         -2.65922977495843870e-001 );
    tau[3] = lapack_make_complex_double( 1.98987975280288490e+000,
                                         -1.41908685396276930e-001 );
}
static void init_c( lapack_int size, lapack_complex_double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    c[0] = lapack_make_complex_double( -3.10981029656006490e-001,
                                       2.62390243772255500e-001 );
    c[8] = lapack_make_complex_double( -3.17534144502360170e-001,
                                       4.83496706343327200e-001 );
    c[16] = lapack_make_complex_double( 4.96614318796456340e-001,
                                        -2.99683439908111190e-001 );
    c[24] = lapack_make_complex_double( -7.19581794453816310e-003,
                                        -3.71789321048053510e-001 );
    c[1] = lapack_make_complex_double( 3.17459801107173260e-001,
                                       -6.41398373665513440e-001 );
    c[9] = lapack_make_complex_double( -2.06186271838591340e-001,
                                       1.57696475525517740e-001 );
    c[17] = lapack_make_complex_double( -7.92590243451013750e-002,
                                        -3.09374948323355990e-001 );
    c[25] = lapack_make_complex_double( -2.81616606005113960e-002,
                                        -1.49146751526478600e-001 );
    c[2] = lapack_make_complex_double( -2.00841914986170850e-001,
                                       1.49011743376836450e-001 );
    c[10] = lapack_make_complex_double( 4.89188100959905730e-001,
                                        -9.00253506243144790e-002 );
    c[18] = lapack_make_complex_double( 3.57457076059659040e-002,
                                        -2.19038212535253160e-002 );
    c[26] = lapack_make_complex_double( 5.62461584914262350e-001,
                                        -7.09935540642336740e-002 );
    c[3] = lapack_make_complex_double( 1.19857271846585860e-001,
                                       -1.23096657572169240e-001 );
    c[11] = lapack_make_complex_double( 2.56601066116824760e-001,
                                        -3.05538478436464810e-001 );
    c[19] = lapack_make_complex_double( 4.48864600443415420e-001,
                                        -2.14082501679258070e-001 );
    c[27] = lapack_make_complex_double( -1.65130169153753950e-001,
                                        1.79976238026742870e-001 );
    c[4] = lapack_make_complex_double( -2.68869015223422270e-001,
                                       -1.65208672004753480e-001 );
    c[12] = lapack_make_complex_double( 1.69670826543481260e-001,
                                        -2.49070210545617840e-001 );
    c[20] = lapack_make_complex_double( -4.95609811295813220e-002,
                                        1.15754472307329690e-001 );
    c[28] = lapack_make_complex_double( -4.88520186437440000e-001,
                                        -4.54037797675943160e-001 );
    c[5] = lapack_make_complex_double( -3.49853658363007300e-001,
                                       9.07028003163352360e-002 );
    c[13] = lapack_make_complex_double( -4.91043340586379550e-002,
                                        -3.13335633697416240e-001 );
    c[21] = lapack_make_complex_double( -1.25647898922387340e-001,
                                        -5.29960507352611470e-001 );
    c[29] = lapack_make_complex_double( 1.03906587226856650e-001,
                                        4.49630621031451990e-002 );
}
static void init_work( lapack_int size, lapack_complex_double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
}

/* Auxiliary function: C interface to zunmbr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zunmbr( lapack_complex_double *c, lapack_complex_double *c_i,
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
