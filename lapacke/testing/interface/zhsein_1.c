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
* zhsein_1 is the test program for the C interface to LAPACK
* routine zhsein
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

static void init_scalars_zhsein( char *job, char *eigsrc, char *initv,
                                 lapack_int *n, lapack_int *ldh,
                                 lapack_int *ldvl, lapack_int *ldvr,
                                 lapack_int *mm );
static void init_select( lapack_int size, lapack_int *select );
static void init_h( lapack_int size, lapack_complex_double *h );
static void init_w( lapack_int size, lapack_complex_double *w );
static void init_vl( lapack_int size, lapack_complex_double *vl );
static void init_vr( lapack_int size, lapack_complex_double *vr );
static void init_work( lapack_int size, lapack_complex_double *work );
static void init_rwork( lapack_int size, double *rwork );
static void init_ifaill( lapack_int size, lapack_int *ifaill );
static void init_ifailr( lapack_int size, lapack_int *ifailr );
static int compare_zhsein( lapack_complex_double *w, lapack_complex_double *w_i,
                           lapack_complex_double *vl,
                           lapack_complex_double *vl_i,
                           lapack_complex_double *vr,
                           lapack_complex_double *vr_i, lapack_int m,
                           lapack_int m_i, lapack_int *ifaill,
                           lapack_int *ifaill_i, lapack_int *ifailr,
                           lapack_int *ifailr_i, lapack_int info,
                           lapack_int info_i, char job, lapack_int ldvl,
                           lapack_int ldvr, lapack_int mm, lapack_int n );

int main(void)
{
    /* Local scalars */
    char job, job_i;
    char eigsrc, eigsrc_i;
    char initv, initv_i;
    lapack_int n, n_i;
    lapack_int ldh, ldh_i;
    lapack_int ldh_r;
    lapack_int ldvl, ldvl_i;
    lapack_int ldvl_r;
    lapack_int ldvr, ldvr_i;
    lapack_int ldvr_r;
    lapack_int mm, mm_i;
    lapack_int m, m_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_int *select = NULL, *select_i = NULL;
    lapack_complex_double *h = NULL, *h_i = NULL;
    lapack_complex_double *w = NULL, *w_i = NULL;
    lapack_complex_double *vl = NULL, *vl_i = NULL;
    lapack_complex_double *vr = NULL, *vr_i = NULL;
    lapack_complex_double *work = NULL, *work_i = NULL;
    double *rwork = NULL, *rwork_i = NULL;
    lapack_int *ifaill = NULL, *ifaill_i = NULL;
    lapack_int *ifailr = NULL, *ifailr_i = NULL;
    lapack_complex_double *w_save = NULL;
    lapack_complex_double *vl_save = NULL;
    lapack_complex_double *vr_save = NULL;
    lapack_int *ifaill_save = NULL;
    lapack_int *ifailr_save = NULL;
    lapack_complex_double *h_r = NULL;
    lapack_complex_double *vl_r = NULL;
    lapack_complex_double *vr_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_zhsein( &job, &eigsrc, &initv, &n, &ldh, &ldvl, &ldvr, &mm );
    ldh_r = n+2;
    ldvl_r = mm+2;
    ldvr_r = mm+2;
    job_i = job;
    eigsrc_i = eigsrc;
    initv_i = initv;
    n_i = n;
    ldh_i = ldh;
    ldvl_i = ldvl;
    ldvr_i = ldvr;
    mm_i = mm;

    /* Allocate memory for the LAPACK routine arrays */
    select = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    h = (lapack_complex_double *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_double) );
    w = (lapack_complex_double *)
        LAPACKE_malloc( n * sizeof(lapack_complex_double) );
    vl = (lapack_complex_double *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_double) );
    vr = (lapack_complex_double *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_double) );
    work = (lapack_complex_double *)
        LAPACKE_malloc( n*n * sizeof(lapack_complex_double) );
    rwork = (double *)LAPACKE_malloc( n * sizeof(double) );
    ifaill = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    h_i = (lapack_complex_double *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_double) );
    w_i = (lapack_complex_double *)
        LAPACKE_malloc( n * sizeof(lapack_complex_double) );
    vl_i = (lapack_complex_double *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_double) );
    vr_i = (lapack_complex_double *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_double) );
    work_i = (lapack_complex_double *)
        LAPACKE_malloc( n*n * sizeof(lapack_complex_double) );
    rwork_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    ifaill_i = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr_i = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    w_save = (lapack_complex_double *)
        LAPACKE_malloc( n * sizeof(lapack_complex_double) );
    vl_save = (lapack_complex_double *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_double) );
    vr_save = (lapack_complex_double *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_double) );
    ifaill_save = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr_save = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the row-major arrays */
    h_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );
    vl_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(mm+2) * sizeof(lapack_complex_double) );
    vr_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(mm+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_select( n, select );
    init_h( ldh*n, h );
    init_w( n, w );
    init_vl( ldvl*mm, vl );
    init_vr( ldvr*mm, vr );
    init_work( n*n, work );
    init_rwork( n, rwork );
    init_ifaill( mm, ifaill );
    init_ifailr( mm, ifailr );

    /* Backup the ouptut arrays */
    for( i = 0; i < n; i++ ) {
        w_save[i] = w[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_save[i] = vl[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_save[i] = vr[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_save[i] = ifaill[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_save[i] = ifailr[i];
    }

    /* Call the LAPACK routine */
    zhsein_( &job, &eigsrc, &initv, select, &n, h, &ldh, w, vl, &ldvl, vr,
             &ldvr, &mm, &m, work, rwork, ifaill, ifailr, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl_save[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr_save[i];
    }
    for( i = 0; i < n*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_i[i] = ifaill_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_i[i] = ifailr_save[i];
    }
    info_i = LAPACKE_zhsein_work( LAPACK_COL_MAJOR, job_i, eigsrc_i, initv_i,
                                  select_i, n_i, h_i, ldh_i, w_i, vl_i, ldvl_i,
                                  vr_i, ldvr_i, mm_i, &m_i, work_i, rwork_i,
                                  ifaill_i, ifailr_i );

    failed = compare_zhsein( w, w_i, vl, vl_i, vr, vr_i, m, m_i, ifaill,
                             ifaill_i, ifailr, ifailr_i, info, info_i, job,
                             ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zhsein\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zhsein\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl_save[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr_save[i];
    }
    for( i = 0; i < n*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_i[i] = ifaill_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_i[i] = ifailr_save[i];
    }
    info_i = LAPACKE_zhsein( LAPACK_COL_MAJOR, job_i, eigsrc_i, initv_i,
                             select_i, n_i, h_i, ldh_i, w_i, vl_i, ldvl_i, vr_i,
                             ldvr_i, mm_i, &m_i, ifaill_i, ifailr_i );

    failed = compare_zhsein( w, w_i, vl, vl_i, vr, vr_i, m, m_i, ifaill,
                             ifaill_i, ifailr, ifailr_i, info, info_i, job,
                             ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zhsein\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zhsein\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl_save[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr_save[i];
    }
    for( i = 0; i < n*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_i[i] = ifaill_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_i[i] = ifailr_save[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_zhsein_work( LAPACK_ROW_MAJOR, job_i, eigsrc_i, initv_i,
                                  select_i, n_i, h_r, ldh_r, w_i, vl_r, ldvl_r,
                                  vr_r, ldvr_r, mm_i, &m_i, work_i, rwork_i,
                                  ifaill_i, ifailr_i );

    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, mm, vl_r, mm+2, vl_i, ldvl );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, mm, vr_r, mm+2, vr_i, ldvr );
    }

    failed = compare_zhsein( w, w_i, vl, vl_i, vr, vr_i, m, m_i, ifaill,
                             ifaill_i, ifailr, ifailr_i, info, info_i, job,
                             ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zhsein\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zhsein\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl_save[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr_save[i];
    }
    for( i = 0; i < n*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_i[i] = ifaill_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_i[i] = ifailr_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_zhsein( LAPACK_ROW_MAJOR, job_i, eigsrc_i, initv_i,
                             select_i, n_i, h_r, ldh_r, w_i, vl_r, ldvl_r, vr_r,
                             ldvr_r, mm_i, &m_i, ifaill_i, ifailr_i );

    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, mm, vl_r, mm+2, vl_i, ldvl );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, mm, vr_r, mm+2, vr_i, ldvr );
    }

    failed = compare_zhsein( w, w_i, vl, vl_i, vr, vr_i, m, m_i, ifaill,
                             ifaill_i, ifailr, ifailr_i, info, info_i, job,
                             ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zhsein\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zhsein\n" );
    }

    /* Release memory */
    if( select != NULL ) {
        LAPACKE_free( select );
    }
    if( select_i != NULL ) {
        LAPACKE_free( select_i );
    }
    if( h != NULL ) {
        LAPACKE_free( h );
    }
    if( h_i != NULL ) {
        LAPACKE_free( h_i );
    }
    if( h_r != NULL ) {
        LAPACKE_free( h_r );
    }
    if( w != NULL ) {
        LAPACKE_free( w );
    }
    if( w_i != NULL ) {
        LAPACKE_free( w_i );
    }
    if( w_save != NULL ) {
        LAPACKE_free( w_save );
    }
    if( vl != NULL ) {
        LAPACKE_free( vl );
    }
    if( vl_i != NULL ) {
        LAPACKE_free( vl_i );
    }
    if( vl_r != NULL ) {
        LAPACKE_free( vl_r );
    }
    if( vl_save != NULL ) {
        LAPACKE_free( vl_save );
    }
    if( vr != NULL ) {
        LAPACKE_free( vr );
    }
    if( vr_i != NULL ) {
        LAPACKE_free( vr_i );
    }
    if( vr_r != NULL ) {
        LAPACKE_free( vr_r );
    }
    if( vr_save != NULL ) {
        LAPACKE_free( vr_save );
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
    if( ifaill != NULL ) {
        LAPACKE_free( ifaill );
    }
    if( ifaill_i != NULL ) {
        LAPACKE_free( ifaill_i );
    }
    if( ifaill_save != NULL ) {
        LAPACKE_free( ifaill_save );
    }
    if( ifailr != NULL ) {
        LAPACKE_free( ifailr );
    }
    if( ifailr_i != NULL ) {
        LAPACKE_free( ifailr_i );
    }
    if( ifailr_save != NULL ) {
        LAPACKE_free( ifailr_save );
    }

    return 0;
}

/* Auxiliary function: zhsein scalar parameters initialization */
static void init_scalars_zhsein( char *job, char *eigsrc, char *initv,
                                 lapack_int *n, lapack_int *ldh,
                                 lapack_int *ldvl, lapack_int *ldvr,
                                 lapack_int *mm )
{
    *job = 'R';
    *eigsrc = 'Q';
    *initv = 'N';
    *n = 4;
    *ldh = 8;
    *ldvl = 8;
    *ldvr = 8;
    *mm = 4;

    return;
}

/* Auxiliary functions: zhsein array parameters initialization */
static void init_select( lapack_int size, lapack_int *select ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        select[i] = 0;
    }
    select[0] = -1;
    select[1] = -1;
    select[2] = 0;
    select[3] = 0;
}
static void init_h( lapack_int size, lapack_complex_double *h ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        h[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    h[0] = lapack_make_complex_double( -3.97000000000000020e+000,
                                       -5.04000000000000000e+000 );
    h[8] = lapack_make_complex_double( -1.13180518733977030e+000,
                                       -2.56930489882743900e+000 );
    h[16] = lapack_make_complex_double( -4.60274243753355350e+000,
                                        -1.42631904083292180e-001 );
    h[24] = lapack_make_complex_double( -1.42491228936652710e+000,
                                        1.73298370334218620e+000 );
    h[1] = lapack_make_complex_double( -5.47965327370263560e+000,
                                       0.00000000000000000e+000 );
    h[9] = lapack_make_complex_double( 1.85847282076558700e+000,
                                       -1.55018070644028950e+000 );
    h[17] = lapack_make_complex_double( 4.41446552691701300e+000,
                                        -7.63823711555098320e-001 );
    h[25] = lapack_make_complex_double( -4.80526133699015420e-001,
                                        -1.19759999733274710e+000 );
    h[2] = lapack_make_complex_double( 6.93222211814628180e-001,
                                       -4.82875276260254950e-001 );
    h[10] = lapack_make_complex_double( 6.26727681806422240e+000,
                                        0.00000000000000000e+000 );
    h[18] = lapack_make_complex_double( -4.50380940334500930e-001,
                                        -2.89818325981801020e-002 );
    h[26] = lapack_make_complex_double( -1.34668445007873290e+000,
                                        1.65792489538873020e+000 );
    h[3] = lapack_make_complex_double( -2.11294690792069330e-001,
                                       8.64412259893682090e-002 );
    h[11] = lapack_make_complex_double( 1.24214618876649560e-001,
                                        -2.28927604979682810e-001 );
    h[19] = lapack_make_complex_double( -3.49998583739325890e+000,
                                        0.00000000000000000e+000 );
    h[27] = lapack_make_complex_double( 2.56190811956891370e+000,
                                        -3.37083746096152880e+000 );
}
static void init_w( lapack_int size, lapack_complex_double *w ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        w[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    w[0] = lapack_make_complex_double( -6.00042534294925110e+000,
                                       -6.99984337157039070e+000 );
    w[1] = lapack_make_complex_double( -5.00003345759696490e+000,
                                       2.00602716231651220e+000 );
    w[2] = lapack_make_complex_double( 7.99819451620824480e+000,
                                       -9.96365091392899080e-001 );
    w[3] = lapack_make_complex_double( 3.00226428433797170e+000,
                                       -3.99981869935322360e+000 );
}
static void init_vl( lapack_int size, lapack_complex_double *vl ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vl[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    vl[0] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[8] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[16] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[24] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[9] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[17] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[25] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[10] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[18] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[26] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[11] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[19] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[27] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
}
static void init_vr( lapack_int size, lapack_complex_double *vr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vr[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    vr[0] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[8] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[16] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[24] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[9] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[17] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[25] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[10] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[18] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[26] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[11] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[19] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[27] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
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
static void init_ifaill( lapack_int size, lapack_int *ifaill ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ifaill[i] = 0;
    }
}
static void init_ifailr( lapack_int size, lapack_int *ifailr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ifailr[i] = 0;
    }
}

/* Auxiliary function: C interface to zhsein results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zhsein( lapack_complex_double *w, lapack_complex_double *w_i,
                           lapack_complex_double *vl,
                           lapack_complex_double *vl_i,
                           lapack_complex_double *vr,
                           lapack_complex_double *vr_i, lapack_int m,
                           lapack_int m_i, lapack_int *ifaill,
                           lapack_int *ifaill_i, lapack_int *ifailr,
                           lapack_int *ifailr_i, lapack_int info,
                           lapack_int info_i, char job, lapack_int ldvl,
                           lapack_int ldvr, lapack_int mm, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < n; i++ ) {
        failed += compare_complex_doubles(w[i],w_i[i]);
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        for( i = 0; i < ldvl*mm; i++ ) {
            failed += compare_complex_doubles(vl[i],vl_i[i]);
        }
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        for( i = 0; i < ldvr*mm; i++ ) {
            failed += compare_complex_doubles(vr[i],vr_i[i]);
        }
    }
    failed += (m == m_i) ? 0 : 1;
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        for( i = 0; i < mm; i++ ) {
            failed += (ifaill[i] == ifaill_i[i]) ? 0 : 1;
        }
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        for( i = 0; i < mm; i++ ) {
            failed += (ifailr[i] == ifailr_i[i]) ? 0 : 1;
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
