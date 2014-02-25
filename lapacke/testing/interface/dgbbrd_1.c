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
* dgbbrd_1 is the test program for the C interface to LAPACK
* routine dgbbrd
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

static void init_scalars_dgbbrd( char *vect, lapack_int *m, lapack_int *n,
                                 lapack_int *ncc, lapack_int *kl,
                                 lapack_int *ku, lapack_int *ldab,
                                 lapack_int *ldq, lapack_int *ldpt,
                                 lapack_int *ldc );
static void init_ab( lapack_int size, double *ab );
static void init_d( lapack_int size, double *d );
static void init_e( lapack_int size, double *e );
static void init_q( lapack_int size, double *q );
static void init_pt( lapack_int size, double *pt );
static void init_c( lapack_int size, double *c );
static void init_work( lapack_int size, double *work );
static int compare_dgbbrd( double *ab, double *ab_i, double *d, double *d_i,
                           double *e, double *e_i, double *q, double *q_i,
                           double *pt, double *pt_i, double *c, double *c_i,
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
    double *ab = NULL, *ab_i = NULL;
    double *d = NULL, *d_i = NULL;
    double *e = NULL, *e_i = NULL;
    double *q = NULL, *q_i = NULL;
    double *pt = NULL, *pt_i = NULL;
    double *c = NULL, *c_i = NULL;
    double *work = NULL, *work_i = NULL;
    double *ab_save = NULL;
    double *d_save = NULL;
    double *e_save = NULL;
    double *q_save = NULL;
    double *pt_save = NULL;
    double *c_save = NULL;
    double *ab_r = NULL;
    double *q_r = NULL;
    double *pt_r = NULL;
    double *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dgbbrd( &vect, &m, &n, &ncc, &kl, &ku, &ldab, &ldq, &ldpt,
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
    ab = (double *)LAPACKE_malloc( ldab*n * sizeof(double) );
    d = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );
    e = (double *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(double) );
    q = (double *)LAPACKE_malloc( ldq*m * sizeof(double) );
    pt = (double *)LAPACKE_malloc( ldpt*n * sizeof(double) );
    c = (double *)LAPACKE_malloc( ldc*ncc * sizeof(double) );
    work = (double *)LAPACKE_malloc( ((2*MAX(m,n))) * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (double *)LAPACKE_malloc( ldab*n * sizeof(double) );
    d_i = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );
    e_i = (double *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(double) );
    q_i = (double *)LAPACKE_malloc( ldq*m * sizeof(double) );
    pt_i = (double *)LAPACKE_malloc( ldpt*n * sizeof(double) );
    c_i = (double *)LAPACKE_malloc( ldc*ncc * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( ((2*MAX(m,n))) * sizeof(double) );

    /* Allocate memory for the backup arrays */
    ab_save = (double *)LAPACKE_malloc( ldab*n * sizeof(double) );
    d_save = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );
    e_save = (double *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(double) );
    q_save = (double *)LAPACKE_malloc( ldq*m * sizeof(double) );
    pt_save = (double *)LAPACKE_malloc( ldpt*n * sizeof(double) );
    c_save = (double *)LAPACKE_malloc( ldc*ncc * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    ab_r = (double *)LAPACKE_malloc( (kl+ku+1)*(n+2) * sizeof(double) );
    q_r = (double *)LAPACKE_malloc( m*(m+2) * sizeof(double) );
    pt_r = (double *)LAPACKE_malloc( n*(n+2) * sizeof(double) );
    c_r = (double *)LAPACKE_malloc( m*(ncc+2) * sizeof(double) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_d( (MIN(m,n)), d );
    init_e( (MIN(m,n)-1), e );
    init_q( ldq*m, q );
    init_pt( ldpt*n, pt );
    init_c( ldc*ncc, c );
    init_work( (2*MAX(m,n)), work );

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
    dgbbrd_( &vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt,
             c, &ldc, work, &info );

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
    for( i = 0; i < (2*MAX(m,n)); i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dgbbrd_work( LAPACK_COL_MAJOR, vect_i, m_i, n_i, ncc_i,
                                  kl_i, ku_i, ab_i, ldab_i, d_i, e_i, q_i,
                                  ldq_i, pt_i, ldpt_i, c_i, ldc_i, work_i );

    failed = compare_dgbbrd( ab, ab_i, d, d_i, e, e_i, q, q_i, pt, pt_i, c, c_i,
                             info, info_i, ldab, ldc, ldpt, ldq, m, n, ncc,
                             vect );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dgbbrd\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dgbbrd\n" );
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
    for( i = 0; i < (2*MAX(m,n)); i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dgbbrd( LAPACK_COL_MAJOR, vect_i, m_i, n_i, ncc_i, kl_i,
                             ku_i, ab_i, ldab_i, d_i, e_i, q_i, ldq_i, pt_i,
                             ldpt_i, c_i, ldc_i );

    failed = compare_dgbbrd( ab, ab_i, d, d_i, e, e_i, q, q_i, pt, pt_i, c, c_i,
                             info, info_i, ldab, ldc, ldpt, ldq, m, n, ncc,
                             vect );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dgbbrd\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dgbbrd\n" );
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
    for( i = 0; i < (2*MAX(m,n)); i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, m, q_i, ldq, q_r, m+2 );
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, pt_i, ldpt, pt_r, n+2 );
    }
    if( ncc != 0 ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, ncc, c_i, ldc, c_r, ncc+2 );
    }
    info_i = LAPACKE_dgbbrd_work( LAPACK_ROW_MAJOR, vect_i, m_i, n_i, ncc_i,
                                  kl_i, ku_i, ab_r, ldab_r, d_i, e_i, q_r,
                                  ldq_r, pt_r, ldpt_r, c_r, ldc_r, work_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, kl+ku+1, n, ab_r, n+2, ab_i, ldab );
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, m, q_r, m+2, q_i, ldq );
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, n, pt_r, n+2, pt_i, ldpt );
    }
    if( ncc != 0 ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, ncc, c_r, ncc+2, c_i, ldc );
    }

    failed = compare_dgbbrd( ab, ab_i, d, d_i, e, e_i, q, q_i, pt, pt_i, c, c_i,
                             info, info_i, ldab, ldc, ldpt, ldq, m, n, ncc,
                             vect );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dgbbrd\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dgbbrd\n" );
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
    for( i = 0; i < (2*MAX(m,n)); i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, m, q_i, ldq, q_r, m+2 );
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, pt_i, ldpt, pt_r, n+2 );
    }
    if( ncc != 0 ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, ncc, c_i, ldc, c_r, ncc+2 );
    }
    info_i = LAPACKE_dgbbrd( LAPACK_ROW_MAJOR, vect_i, m_i, n_i, ncc_i, kl_i,
                             ku_i, ab_r, ldab_r, d_i, e_i, q_r, ldq_r, pt_r,
                             ldpt_r, c_r, ldc_r );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, kl+ku+1, n, ab_r, n+2, ab_i, ldab );
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, m, q_r, m+2, q_i, ldq );
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, n, pt_r, n+2, pt_i, ldpt );
    }
    if( ncc != 0 ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, ncc, c_r, ncc+2, c_i, ldc );
    }

    failed = compare_dgbbrd( ab, ab_i, d, d_i, e, e_i, q, q_i, pt, pt_i, c, c_i,
                             info, info_i, ldab, ldc, ldpt, ldq, m, n, ncc,
                             vect );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dgbbrd\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dgbbrd\n" );
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

    return 0;
}

/* Auxiliary function: dgbbrd scalar parameters initialization */
static void init_scalars_dgbbrd( char *vect, lapack_int *m, lapack_int *n,
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

/* Auxiliary functions: dgbbrd array parameters initialization */
static void init_ab( lapack_int size, double *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = 0;
    }
    ab[0] = 0.00000000000000000e+000;  /* ab[0,0] */
    ab[17] = -1.28000000000000000e+000;  /* ab[0,1] */
    ab[34] = -3.10000000000000000e-001;  /* ab[0,2] */
    ab[51] = -3.49999999999999980e-001;  /* ab[0,3] */
    ab[1] = -5.69999999999999950e-001;  /* ab[1,0] */
    ab[18] = 1.08000000000000010e+000;  /* ab[1,1] */
    ab[35] = 4.00000000000000020e-001;  /* ab[1,2] */
    ab[52] = 8.00000000000000020e-002;  /* ab[1,3] */
    ab[2] = -1.92999999999999990e+000;  /* ab[2,0] */
    ab[19] = 2.39999999999999990e-001;  /* ab[2,1] */
    ab[36] = -6.60000000000000030e-001;  /* ab[2,2] */
    ab[53] = -2.12999999999999990e+000;  /* ab[2,3] */
    ab[3] = 2.29999999999999980e+000;  /* ab[3,0] */
    ab[20] = 6.40000000000000010e-001;  /* ab[3,1] */
    ab[37] = 1.49999999999999990e-001;  /* ab[3,2] */
    ab[54] = 5.00000000000000000e-001;  /* ab[3,3] */
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
static void init_q( lapack_int size, double *q ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        q[i] = 0;
    }
}
static void init_pt( lapack_int size, double *pt ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        pt[i] = 0;
    }
}
static void init_c( lapack_int size, double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = 0;
    }
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to dgbbrd results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dgbbrd( double *ab, double *ab_i, double *d, double *d_i,
                           double *e, double *e_i, double *q, double *q_i,
                           double *pt, double *pt_i, double *c, double *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldab,
                           lapack_int ldc, lapack_int ldpt, lapack_int ldq,
                           lapack_int m, lapack_int n, lapack_int ncc,
                           char vect )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldab*n; i++ ) {
        failed += compare_doubles(ab[i],ab_i[i]);
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        failed += compare_doubles(d[i],d_i[i]);
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        failed += compare_doubles(e[i],e_i[i]);
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'q' ) ) {
        for( i = 0; i < ldq*m; i++ ) {
            failed += compare_doubles(q[i],q_i[i]);
        }
    }
    if( LAPACKE_lsame( vect, 'b' ) || LAPACKE_lsame( vect, 'p' ) ) {
        for( i = 0; i < ldpt*n; i++ ) {
            failed += compare_doubles(pt[i],pt_i[i]);
        }
    }
    if( ncc != 0 ) {
        for( i = 0; i < ldc*ncc; i++ ) {
            failed += compare_doubles(c[i],c_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
