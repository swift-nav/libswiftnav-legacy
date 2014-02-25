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
* chsein_1 is the test program for the C interface to LAPACK
* routine chsein
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

static void init_scalars_chsein( char *job, char *eigsrc, char *initv,
                                 lapack_int *n, lapack_int *ldh,
                                 lapack_int *ldvl, lapack_int *ldvr,
                                 lapack_int *mm );
static void init_select( lapack_int size, lapack_int *select );
static void init_h( lapack_int size, lapack_complex_float *h );
static void init_w( lapack_int size, lapack_complex_float *w );
static void init_vl( lapack_int size, lapack_complex_float *vl );
static void init_vr( lapack_int size, lapack_complex_float *vr );
static void init_work( lapack_int size, lapack_complex_float *work );
static void init_rwork( lapack_int size, float *rwork );
static void init_ifaill( lapack_int size, lapack_int *ifaill );
static void init_ifailr( lapack_int size, lapack_int *ifailr );
static int compare_chsein( lapack_complex_float *w, lapack_complex_float *w_i,
                           lapack_complex_float *vl, lapack_complex_float *vl_i,
                           lapack_complex_float *vr, lapack_complex_float *vr_i,
                           lapack_int m, lapack_int m_i, lapack_int *ifaill,
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
    lapack_complex_float *h = NULL, *h_i = NULL;
    lapack_complex_float *w = NULL, *w_i = NULL;
    lapack_complex_float *vl = NULL, *vl_i = NULL;
    lapack_complex_float *vr = NULL, *vr_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    float *rwork = NULL, *rwork_i = NULL;
    lapack_int *ifaill = NULL, *ifaill_i = NULL;
    lapack_int *ifailr = NULL, *ifailr_i = NULL;
    lapack_complex_float *w_save = NULL;
    lapack_complex_float *vl_save = NULL;
    lapack_complex_float *vr_save = NULL;
    lapack_int *ifaill_save = NULL;
    lapack_int *ifailr_save = NULL;
    lapack_complex_float *h_r = NULL;
    lapack_complex_float *vl_r = NULL;
    lapack_complex_float *vr_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_chsein( &job, &eigsrc, &initv, &n, &ldh, &ldvl, &ldvr, &mm );
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
    h = (lapack_complex_float *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_float) );
    w = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );
    vl = (lapack_complex_float *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_float) );
    vr = (lapack_complex_float *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( n*n * sizeof(lapack_complex_float) );
    rwork = (float *)LAPACKE_malloc( n * sizeof(float) );
    ifaill = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    h_i = (lapack_complex_float *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_float) );
    w_i = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );
    vl_i = (lapack_complex_float *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_float) );
    vr_i = (lapack_complex_float *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( n*n * sizeof(lapack_complex_float) );
    rwork_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    ifaill_i = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr_i = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    w_save = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );
    vl_save = (lapack_complex_float *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_float) );
    vr_save = (lapack_complex_float *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_float) );
    ifaill_save = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr_save = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the row-major arrays */
    h_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );
    vl_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(mm+2) * sizeof(lapack_complex_float) );
    vr_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(mm+2) * sizeof(lapack_complex_float) );

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
    chsein_( &job, &eigsrc, &initv, select, &n, h, &ldh, w, vl, &ldvl, vr,
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
    info_i = LAPACKE_chsein_work( LAPACK_COL_MAJOR, job_i, eigsrc_i, initv_i,
                                  select_i, n_i, h_i, ldh_i, w_i, vl_i, ldvl_i,
                                  vr_i, ldvr_i, mm_i, &m_i, work_i, rwork_i,
                                  ifaill_i, ifailr_i );

    failed = compare_chsein( w, w_i, vl, vl_i, vr, vr_i, m, m_i, ifaill,
                             ifaill_i, ifailr, ifailr_i, info, info_i, job,
                             ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to chsein\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to chsein\n" );
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
    info_i = LAPACKE_chsein( LAPACK_COL_MAJOR, job_i, eigsrc_i, initv_i,
                             select_i, n_i, h_i, ldh_i, w_i, vl_i, ldvl_i, vr_i,
                             ldvr_i, mm_i, &m_i, ifaill_i, ifailr_i );

    failed = compare_chsein( w, w_i, vl, vl_i, vr, vr_i, m, m_i, ifaill,
                             ifaill_i, ifailr, ifailr_i, info, info_i, job,
                             ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to chsein\n" );
    } else {
        printf( "FAILED: column-major high-level interface to chsein\n" );
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

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_chsein_work( LAPACK_ROW_MAJOR, job_i, eigsrc_i, initv_i,
                                  select_i, n_i, h_r, ldh_r, w_i, vl_r, ldvl_r,
                                  vr_r, ldvr_r, mm_i, &m_i, work_i, rwork_i,
                                  ifaill_i, ifailr_i );

    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, mm, vl_r, mm+2, vl_i, ldvl );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, mm, vr_r, mm+2, vr_i, ldvr );
    }

    failed = compare_chsein( w, w_i, vl, vl_i, vr, vr_i, m, m_i, ifaill,
                             ifaill_i, ifailr, ifailr_i, info, info_i, job,
                             ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to chsein\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to chsein\n" );
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
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_chsein( LAPACK_ROW_MAJOR, job_i, eigsrc_i, initv_i,
                             select_i, n_i, h_r, ldh_r, w_i, vl_r, ldvl_r, vr_r,
                             ldvr_r, mm_i, &m_i, ifaill_i, ifailr_i );

    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, mm, vl_r, mm+2, vl_i, ldvl );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, mm, vr_r, mm+2, vr_i, ldvr );
    }

    failed = compare_chsein( w, w_i, vl, vl_i, vr, vr_i, m, m_i, ifaill,
                             ifaill_i, ifailr, ifailr_i, info, info_i, job,
                             ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to chsein\n" );
    } else {
        printf( "FAILED: row-major high-level interface to chsein\n" );
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

/* Auxiliary function: chsein scalar parameters initialization */
static void init_scalars_chsein( char *job, char *eigsrc, char *initv,
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

/* Auxiliary functions: chsein array parameters initialization */
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
static void init_h( lapack_int size, lapack_complex_float *h ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        h[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    h[0] = lapack_make_complex_float( -3.970000029e+000, -5.039999962e+000 );
    h[8] = lapack_make_complex_float( -1.131804943e+000, -2.569304705e+000 );
    h[16] = lapack_make_complex_float( -4.602741241e+000, -1.426316500e-001 );
    h[24] = lapack_make_complex_float( -1.424912333e+000, 1.732983828e+000 );
    h[1] = lapack_make_complex_float( -5.479653358e+000, 0.000000000e+000 );
    h[9] = lapack_make_complex_float( 1.858472466e+000, -1.550180435e+000 );
    h[17] = lapack_make_complex_float( 4.414464474e+000, -7.638237476e-001 );
    h[25] = lapack_make_complex_float( -4.805260897e-001, -1.197600603e+000 );
    h[2] = lapack_make_complex_float( 6.932221651e-001, -4.828752279e-001 );
    h[10] = lapack_make_complex_float( 6.267275810e+000, 0.000000000e+000 );
    h[18] = lapack_make_complex_float( -4.503800869e-001, -2.898204327e-002 );
    h[26] = lapack_make_complex_float( -1.346683741e+000, 1.657924891e+000 );
    h[3] = lapack_make_complex_float( -2.112946808e-001, 8.644121885e-002 );
    h[11] = lapack_make_complex_float( 1.242147088e-001, -2.289276123e-001 );
    h[19] = lapack_make_complex_float( -3.499985933e+000, 0.000000000e+000 );
    h[27] = lapack_make_complex_float( 2.561908484e+000, -3.370837450e+000 );
}
static void init_w( lapack_int size, lapack_complex_float *w ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        w[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    w[0] = lapack_make_complex_float( -6.000423908e+000, -6.999841690e+000 );
    w[1] = lapack_make_complex_float( -5.000031471e+000, 2.006026745e+000 );
    w[2] = lapack_make_complex_float( 7.998193264e+000, -9.963648915e-001 );
    w[3] = lapack_make_complex_float( 3.002264738e+000, -3.999819279e+000 );
}
static void init_vl( lapack_int size, lapack_complex_float *vl ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vl[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    vl[0] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[8] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[16] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[24] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[9] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[17] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[25] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[10] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[18] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[26] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[19] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[27] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_vr( lapack_int size, lapack_complex_float *vr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vr[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    vr[0] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[8] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[16] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[24] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[9] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[17] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[25] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[10] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[18] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[26] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[19] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[27] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}
static void init_rwork( lapack_int size, float *rwork ) {
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

/* Auxiliary function: C interface to chsein results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_chsein( lapack_complex_float *w, lapack_complex_float *w_i,
                           lapack_complex_float *vl, lapack_complex_float *vl_i,
                           lapack_complex_float *vr, lapack_complex_float *vr_i,
                           lapack_int m, lapack_int m_i, lapack_int *ifaill,
                           lapack_int *ifaill_i, lapack_int *ifailr,
                           lapack_int *ifailr_i, lapack_int info,
                           lapack_int info_i, char job, lapack_int ldvl,
                           lapack_int ldvr, lapack_int mm, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < n; i++ ) {
        failed += compare_complex_floats(w[i],w_i[i]);
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        for( i = 0; i < ldvl*mm; i++ ) {
            failed += compare_complex_floats(vl[i],vl_i[i]);
        }
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        for( i = 0; i < ldvr*mm; i++ ) {
            failed += compare_complex_floats(vr[i],vr_i[i]);
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
