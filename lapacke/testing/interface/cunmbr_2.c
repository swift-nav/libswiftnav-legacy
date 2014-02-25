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
* cunmbr_2 is the test program for the C interface to LAPACK
* routine cunmbr
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

static void init_scalars_cunmbr( char *vect, char *side, char *trans,
                                 lapack_int *m, lapack_int *n, lapack_int *k,
                                 lapack_int *lda, lapack_int *ldc,
                                 lapack_int *lwork );
static void init_a( lapack_int size, lapack_complex_float *a );
static void init_tau( lapack_int size, lapack_complex_float *tau );
static void init_c( lapack_int size, lapack_complex_float *c );
static void init_work( lapack_int size, lapack_complex_float *work );
static int compare_cunmbr( lapack_complex_float *c, lapack_complex_float *c_i,
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
    lapack_complex_float *a = NULL, *a_i = NULL;
    lapack_complex_float *tau = NULL, *tau_i = NULL;
    lapack_complex_float *c = NULL, *c_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    lapack_complex_float *c_save = NULL;
    lapack_complex_float *a_r = NULL;
    lapack_complex_float *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_cunmbr( &vect, &side, &trans, &m, &n, &k, &lda, &ldc, &lwork );
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
    a = (lapack_complex_float *)
        LAPACKE_malloc( (lda*(MIN(nq,k))) * sizeof(lapack_complex_float) );
    tau = (lapack_complex_float *)
        LAPACKE_malloc( MIN(nq,k) * sizeof(lapack_complex_float) );
    c = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_float *)
        LAPACKE_malloc( (lda*(MIN(nq,k))) * sizeof(lapack_complex_float) );
    tau_i = (lapack_complex_float *)
        LAPACKE_malloc( MIN(nq,k) * sizeof(lapack_complex_float) );
    c_i = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    c_save = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_float *)
        LAPACKE_malloc( (r*(MIN(nq,k)+2)) * sizeof(lapack_complex_float) );
    c_r = (lapack_complex_float *)
        LAPACKE_malloc( m*(n+2) * sizeof(lapack_complex_float) );

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
    cunmbr_( &vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work,
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
    info_i = LAPACKE_cunmbr_work( LAPACK_COL_MAJOR, vect_i, side_i, trans_i,
                                  m_i, n_i, k_i, a_i, lda_i, tau_i, c_i, ldc_i,
                                  work_i, lwork_i );

    failed = compare_cunmbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to cunmbr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to cunmbr\n" );
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
    info_i = LAPACKE_cunmbr( LAPACK_COL_MAJOR, vect_i, side_i, trans_i, m_i,
                             n_i, k_i, a_i, lda_i, tau_i, c_i, ldc_i );

    failed = compare_cunmbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to cunmbr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to cunmbr\n" );
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

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, r, MIN(nq, k ), a_i, lda, a_r, MIN(nq,
                       k)+2);
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_cunmbr_work( LAPACK_ROW_MAJOR, vect_i, side_i, trans_i,
                                  m_i, n_i, k_i, a_r, lda_r, tau_i, c_r, ldc_r,
                                  work_i, lwork_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_cunmbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to cunmbr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to cunmbr\n" );
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
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, r, MIN(nq, k ), a_i, lda, a_r, MIN(nq,
                       k)+2);
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_cunmbr( LAPACK_ROW_MAJOR, vect_i, side_i, trans_i, m_i,
                             n_i, k_i, a_r, lda_r, tau_i, c_r, ldc_r );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_cunmbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to cunmbr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to cunmbr\n" );
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

/* Auxiliary function: cunmbr scalar parameters initialization */
static void init_scalars_cunmbr( char *vect, char *side, char *trans,
                                 lapack_int *m, lapack_int *n, lapack_int *k,
                                 lapack_int *lda, lapack_int *ldc,
                                 lapack_int *lwork )
{
    *vect = 'P';
    *side = 'L';
    *trans = 'C';
    *m = 3;
    *n = 4;
    *k = 3;
    *lda = 8;
    *ldc = 8;
    *lwork = 1024;

    return;
}

/* Auxiliary functions: cunmbr array parameters initialization */
static void init_a( lapack_int size, lapack_complex_float *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    a[0] = lapack_make_complex_float( 2.761475801e+000, 0.000000000e+000 );
    a[8] = lapack_make_complex_float( -9.499514103e-001, 0.000000000e+000 );
    a[16] = lapack_make_complex_float( 7.651669532e-002, -2.179353982e-001 );
    a[1] = lapack_make_complex_float( -1.645794213e-001, -2.483377308e-001 );
    a[9] = lapack_make_complex_float( 1.629762888e+000, 0.000000000e+000 );
    a[17] = lapack_make_complex_float( -1.018280864e+000, 0.000000000e+000 );
    a[2] = lapack_make_complex_float( -2.072303614e-004, 1.368010789e-001 );
    a[10] = lapack_make_complex_float( 2.649169862e-001, -2.632617950e-001 );
    a[18] = lapack_make_complex_float( -1.327486753e+000, 0.000000000e+000 );
}
static void init_tau( lapack_int size, lapack_complex_float *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    tau[0] = lapack_make_complex_float( 1.688534141e+000, 5.957157016e-001 );
    tau[1] = lapack_make_complex_float( 1.935088515e+000, 3.544141948e-001 );
    tau[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_c( lapack_int size, lapack_complex_float *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    c[0] = lapack_make_complex_float( -1.258138418e-001, 1.617605835e-001 );
    c[8] = lapack_make_complex_float( -2.246674895e-001, 3.864280879e-001 );
    c[16] = lapack_make_complex_float( 3.459879458e-001, 2.156807780e-001 );
    c[24] = lapack_make_complex_float( -7.099492550e-001, -2.965611219e-001 );
    c[1] = lapack_make_complex_float( -1.163461953e-001, -6.379657388e-001 );
    c[9] = lapack_make_complex_float( -3.240494728e-001, 4.271534085e-001 );
    c[17] = lapack_make_complex_float( -1.995496750e-001, -5.008659363e-001 );
    c[25] = lapack_make_complex_float( -3.233429790e-002, -1.620414853e-002 );
    c[2] = lapack_make_complex_float( -4.606565833e-001, 1.090034544e-001 );
    c[10] = lapack_make_complex_float( 2.170982510e-001, -4.061889052e-001 );
    c[18] = lapack_make_complex_float( 2.733075619e-001, -6.106228232e-001 );
    c[26] = lapack_make_complex_float( -9.940487146e-002, -3.261196911e-001 );
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}

/* Auxiliary function: C interface to cunmbr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_cunmbr( lapack_complex_float *c, lapack_complex_float *c_i,
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
