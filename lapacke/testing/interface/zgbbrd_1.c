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
* zgbbrd_1 is the test program for the C interface to LAPACK
* routine zgbbrd
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

static void init_scalars_zgbbrd( char *vect, lapack_int *m, lapack_int *n,
                                 lapack_int *ncc, lapack_int *kl,
                                 lapack_int *ku, lapack_int *ldab,
                                 lapack_int *ldq, lapack_int *ldpt,
                                 lapack_int *ldc );
static void init_ab( lapack_int size, lapack_complex_double *ab );
static void init_d( lapack_int size, double *d );
static void init_e( lapack_int size, double *e );
static void init_q( lapack_int size, lapack_complex_double *q );
static void init_pt( lapack_int size, lapack_complex_double *pt );
static void init_c( lapack_int size, lapack_complex_double *c );
static void init_work( lapack_int size, lapack_complex_double *work );
static void init_rwork( lapack_int size, double *rwork );
static int compare_zgbbrd( lapack_complex_double *ab,
                           lapack_complex_double *ab_i, double *d, double *d_i,
                           double *e, double *e_i, lapack_complex_double *q,
                           lapack_complex_double *q_i,
                           lapack_complex_double *pt,
                           lapack_complex_double *pt_i,
                           lapack_complex_double *c, lapack_complex_double *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldab,
                           lapack_int ldc, lapack_int ldpt, lapack_int ldq,
                           lapack_int m, lapack_int n, lapack_int ncc,
                           char vect );

int main(void)
{
    /* Local scalars */
    char vect, vect_i;
    lapack_int m, m_i;
    lapack_int n, n_i;
    lapack_int ncc, ncc_i;
    lapack_int kl, kl_i;
    lapack_int ku, ku_i;
    lapack_int ldab, ldab_i;
    lapack_int ldab_r;
    lapack_int ldq, ldq_i;
    lapack_int ldq_r;
    lapack_int ldpt, ldpt_i;
    lapack_int ldpt_r;
    lapack_int ldc, ldc_i;
    lapack_int ldc_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_double *ab = NULL, *ab_i = NULL;
    double *d = NULL, *d_i = NULL;
    double *e = NULL, *e_i = NULL;
    lapack_complex_double *q = NULL, *q_i = NULL;
    lapack_complex_double *pt = NULL, *pt_i = NULL;
    lapack_complex_double *c = NULL, *c_i = NULL;
    lapack_complex_double *work = NULL, *work_i = NULL;
    double *rwork = NULL, *rwork_i = NULL;
    lapack_complex_double *ab_save = NULL;
    double *d_save = NULL;
    double *e_save = NULL;
    lapack_complex_double *q_save = NULL;
    lapack_complex_double *pt_save = NULL;
    lapack_complex_double *c_save = NULL;
    lapack_complex_double *ab_r = NULL;
    lapack_complex_double *q_r = NULL;
    lapack_complex_double *pt_r = NULL;
    lapack_complex_double *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_zgbbrd( &vect, &m, &n, &ncc, &kl, &ku, &ldab, &ldq, &ldpt,
                         &ldc );
    ldab_r = n+2;
    ldq_r = m+2;
    ldpt_r = n+2;
    ldc_r = ncc+2;
    vect_i = vect;
    m_i = m;
    n_i = n;
    ncc_i = ncc;
    kl_i = kl;
    ku_i = ku;
    ldab_i = ldab;
    ldq_i = ldq;
    ldpt_i = ldpt;
    ldc_i = ldc;

    /* Allocate memory for the LAPACK routine arrays */
    ab = (lapack_complex_double *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_double) );
    d = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );
    e = (double *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(double) );
    q = (lapack_complex_double *)
        LAPACKE_malloc( ldq*m * sizeof(lapack_complex_double) );
    pt = (lapack_complex_double *)
        LAPACKE_malloc( ldpt*n * sizeof(lapack_complex_double) );
    c = (lapack_complex_double *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_double) );
    work = (lapack_complex_double *)
        LAPACKE_malloc( MAX(m,n) * sizeof(lapack_complex_double) );
    rwork = (double *)LAPACKE_malloc( MAX(m,n) * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (lapack_complex_double *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_double) );
    d_i = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );
    e_i = (double *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(double) );
    q_i = (lapack_complex_double *)
        LAPACKE_malloc( ldq*m * sizeof(lapack_complex_double) );
    pt_i = (lapack_complex_double *)
        LAPACKE_malloc( ldpt*n * sizeof(lapack_complex_double) );
    c_i = (lapack_complex_double *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_double) );
    work_i = (lapack_complex_double *)
        LAPACKE_malloc( MAX(m,n) * sizeof(lapack_complex_double) );
    rwork_i = (double *)LAPACKE_malloc( MAX(m,n) * sizeof(double) );

    /* Allocate memory for the backup arrays */
    ab_save = (lapack_complex_double *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_double) );
    d_save = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );
    e_save = (double *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(double) );
    q_save = (lapack_complex_double *)
        LAPACKE_malloc( ldq*m * sizeof(lapack_complex_double) );
    pt_save = (lapack_complex_double *)
        LAPACKE_malloc( ldpt*n * sizeof(lapack_complex_double) );
    c_save = (lapack_complex_double *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_double) );

    /* Allocate memory for the row-major arrays */
    ab_r = (lapack_complex_double *)
        LAPACKE_malloc( (kl+ku+1)*(n+2) * sizeof(lapack_complex_double) );
    q_r = (lapack_complex_double *)
        LAPACKE_malloc( m*(m+2) * sizeof(lapack_complex_double) );
    pt_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );
    c_r = (lapack_complex_double *)
        LAPACKE_malloc( m*(ncc+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_d( (MIN(m,n)), d );
    init_e( (MIN(m,n)-1), e );
    init_q( ldq*m, q );
    init_pt( ldpt*n, pt );
    init_c( ldc*ncc, c );
    init_work( (MAX(m,n)), work );
    init_rwork( (MAX(m,n)), rwork );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldab*n; i++ ) {
        ab_save[i] = ab[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_save[i] = d[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_save[i] = e[i];
    }
    for( i = 0; i < ldq*m; i++ ) {
        q_save[i] = q[i];
    }
    for( i = 0; i < ldpt*n; i++ ) {
        pt_save[i] = pt[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    zgbbrd_( &vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt,
             c, &ldc, work, rwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldq*m; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < ldpt*n; i++ ) {
        pt_i[i] = pt_save[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < (MAX(m,n)); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < (MAX(m,n)); i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_zgbbrd_work( LAPACK_COL_MAJOR, vect_i, m_i, n_i, ncc_i,
                                  kl_i, ku_i, ab_i, ldab_i, d_i, e_i, q_i,
                                  ldq_i, pt_i, ldpt_i, c_i, ldc_i, work_i,
                                  rwork_i );

    failed = compare_zgbbrd( ab, ab_i, d, d_i, e, e_i, q, q_i, pt, pt_i, c, c_i,
                             info, info_i, ldab, ldc, ldpt, ldq, m, n, ncc,
                             vect );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zgbbrd\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zgbbrd\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldq*m; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < ldpt*n; i++ ) {
        pt_i[i] = pt_save[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < (MAX(m,n)); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < (MAX(m,n)); i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_zgbbrd( LAPACK_COL_MAJOR, vect_i, m_i, n_i, ncc_i, kl_i,
                             ku_i, ab_i, ldab_i, d_i, e_i, q_i, ldq_i, pt_i,
                             ldpt_i, c_i, ldc_i );

    failed = compare_zgbbrd( ab, ab_i, d, d_i, e, e_i, q, q_i, pt, pt_i, c, c_i,
                             info, info_i, ldab, ldc, ldpt, ldq, m, n, ncc,
                             vect );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zgbbrd\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zgbbrd\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldq*m; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < ldpt*n; i++ ) {
        pt_i[i] = pt_save[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < (MAX(m,n)); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < (MAX(m,n)); i++ ) {
        rwork_i[i] = rwork[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, m, q_i, ldq, q_r, m+2 );
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, pt_i, ldpt, pt_r, n+2 );
    }
    if( ncc != 0 ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, ncc, c_i, ldc, c_r, ncc+2 );
    }
    info_i = LAPACKE_zgbbrd_work( LAPACK_ROW_MAJOR, vect_i, m_i, n_i, ncc_i,
                                  kl_i, ku_i, ab_r, ldab_r, d_i, e_i, q_r,
                                  ldq_r, pt_r, ldpt_r, c_r, ldc_r, work_i,
                                  rwork_i );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, kl+ku+1, n, ab_r, n+2, ab_i, ldab );
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, m, q_r, m+2, q_i, ldq );
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, n, pt_r, n+2, pt_i, ldpt );
    }
    if( ncc != 0 ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, ncc, c_r, ncc+2, c_i, ldc );
    }

    failed = compare_zgbbrd( ab, ab_i, d, d_i, e, e_i, q, q_i, pt, pt_i, c, c_i,
                             info, info_i, ldab, ldc, ldpt, ldq, m, n, ncc,
                             vect );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zgbbrd\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zgbbrd\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldq*m; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < ldpt*n; i++ ) {
        pt_i[i] = pt_save[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < (MAX(m,n)); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < (MAX(m,n)); i++ ) {
        rwork_i[i] = rwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, m, q_i, ldq, q_r, m+2 );
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, pt_i, ldpt, pt_r, n+2 );
    }
    if( ncc != 0 ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, ncc, c_i, ldc, c_r, ncc+2 );
    }
    info_i = LAPACKE_zgbbrd( LAPACK_ROW_MAJOR, vect_i, m_i, n_i, ncc_i, kl_i,
                             ku_i, ab_r, ldab_r, d_i, e_i, q_r, ldq_r, pt_r,
                             ldpt_r, c_r, ldc_r );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, kl+ku+1, n, ab_r, n+2, ab_i, ldab );
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, m, q_r, m+2, q_i, ldq );
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, n, pt_r, n+2, pt_i, ldpt );
    }
    if( ncc != 0 ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, ncc, c_r, ncc+2, c_i, ldc );
    }

    failed = compare_zgbbrd( ab, ab_i, d, d_i, e, e_i, q, q_i, pt, pt_i, c, c_i,
                             info, info_i, ldab, ldc, ldpt, ldq, m, n, ncc,
                             vect );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zgbbrd\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zgbbrd\n" );
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
    if( d != NULL ) {
        LAPACKE_free( d );
    }
    if( d_i != NULL ) {
        LAPACKE_free( d_i );
    }
    if( d_save != NULL ) {
        LAPACKE_free( d_save );
    }
    if( e != NULL ) {
        LAPACKE_free( e );
    }
    if( e_i != NULL ) {
        LAPACKE_free( e_i );
    }
    if( e_save != NULL ) {
        LAPACKE_free( e_save );
    }
    if( q != NULL ) {
        LAPACKE_free( q );
    }
    if( q_i != NULL ) {
        LAPACKE_free( q_i );
    }
    if( q_r != NULL ) {
        LAPACKE_free( q_r );
    }
    if( q_save != NULL ) {
        LAPACKE_free( q_save );
    }
    if( pt != NULL ) {
        LAPACKE_free( pt );
    }
    if( pt_i != NULL ) {
        LAPACKE_free( pt_i );
    }
    if( pt_r != NULL ) {
        LAPACKE_free( pt_r );
    }
    if( pt_save != NULL ) {
        LAPACKE_free( pt_save );
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
    if( rwork != NULL ) {
        LAPACKE_free( rwork );
    }
    if( rwork_i != NULL ) {
        LAPACKE_free( rwork_i );
    }

    return 0;
}

/* Auxiliary function: zgbbrd scalar parameters initialization */
static void init_scalars_zgbbrd( char *vect, lapack_int *m, lapack_int *n,
                                 lapack_int *ncc, lapack_int *kl,
                                 lapack_int *ku, lapack_int *ldab,
                                 lapack_int *ldq, lapack_int *ldpt,
                                 lapack_int *ldc )
{
    *vect = 'N';
    *m = 6;
    *n = 4;
    *ncc = 0;
    *kl = 2;
    *ku = 1;
    *ldab = 17;
    *ldq = 8;
    *ldpt = 8;
    *ldc = 8;

    return;
}

/* Auxiliary functions: zgbbrd array parameters initialization */
static void init_ab( lapack_int size, lapack_complex_double *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    ab[0] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    ab[17] = lapack_make_complex_double( -2.99999999999999990e-002,
                                         9.59999999999999960e-001 );
    ab[34] = lapack_make_complex_double( -6.60000000000000030e-001,
                                         4.19999999999999980e-001 );
    ab[51] = lapack_make_complex_double( -1.11000000000000010e+000,
                                         5.99999999999999980e-001 );
    ab[1] = lapack_make_complex_double( 9.59999999999999960e-001,
                                        -8.10000000000000050e-001 );
    ab[18] = lapack_make_complex_double( -1.20000000000000000e+000,
                                         1.90000000000000000e-001 );
    ab[35] = lapack_make_complex_double( 6.30000000000000000e-001,
                                         -1.70000000000000010e-001 );
    ab[52] = lapack_make_complex_double( 2.20000000000000000e-001,
                                         -2.00000000000000010e-001 );
    ab[2] = lapack_make_complex_double( -9.79999999999999980e-001,
                                        1.98000000000000000e+000 );
    ab[19] = lapack_make_complex_double( 1.01000000000000000e+000,
                                         2.00000000000000000e-002 );
    ab[36] = lapack_make_complex_double( -9.79999999999999980e-001,
                                         -3.59999999999999990e-001 );
    ab[53] = lapack_make_complex_double( 1.47000000000000000e+000,
                                         1.59000000000000010e+000 );
    ab[3] = lapack_make_complex_double( 6.20000000000000000e-001,
                                        -4.60000000000000020e-001 );
    ab[20] = lapack_make_complex_double( 1.90000000000000000e-001,
                                         -5.40000000000000040e-001 );
    ab[37] = lapack_make_complex_double( -1.70000000000000010e-001,
                                         -4.60000000000000020e-001 );
    ab[54] = lapack_make_complex_double( 2.60000000000000010e-001,
                                         2.60000000000000010e-001 );
}
static void init_d( lapack_int size, double *d ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        d[i] = 0;
    }
}
static void init_e( lapack_int size, double *e ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        e[i] = 0;
    }
}
static void init_q( lapack_int size, lapack_complex_double *q ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        q[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
}
static void init_pt( lapack_int size, lapack_complex_double *pt ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        pt[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
}
static void init_c( lapack_int size, lapack_complex_double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
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

/* Auxiliary function: C interface to zgbbrd results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zgbbrd( lapack_complex_double *ab,
                           lapack_complex_double *ab_i, double *d, double *d_i,
                           double *e, double *e_i, lapack_complex_double *q,
                           lapack_complex_double *q_i,
                           lapack_complex_double *pt,
                           lapack_complex_double *pt_i,
                           lapack_complex_double *c, lapack_complex_double *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldab,
                           lapack_int ldc, lapack_int ldpt, lapack_int ldq,
                           lapack_int m, lapack_int n, lapack_int ncc,
                           char vect )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldab*n; i++ ) {
        failed += compare_complex_doubles(ab[i],ab_i[i]);
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        failed += compare_doubles(d[i],d_i[i]);
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        failed += compare_doubles(e[i],e_i[i]);
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        for( i = 0; i < ldq*m; i++ ) {
            failed += compare_complex_doubles(q[i],q_i[i]);
        }
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        for( i = 0; i < ldpt*n; i++ ) {
            failed += compare_complex_doubles(pt[i],pt_i[i]);
        }
    }
    if( ncc != 0 ) {
        for( i = 0; i < ldc*ncc; i++ ) {
            failed += compare_complex_doubles(c[i],c_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
