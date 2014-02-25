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
* ztrsna_1 is the test program for the C interface to LAPACK
* routine ztrsna
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

static void init_scalars_ztrsna( char *job, char *howmny, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldvl,
                                 lapack_int *ldvr, lapack_int *mm,
                                 lapack_int *ldwork );
static void init_select( lapack_int size, lapack_int *select );
static void init_t( lapack_int size, lapack_complex_double *t );
static void init_vl( lapack_int size, lapack_complex_double *vl );
static void init_vr( lapack_int size, lapack_complex_double *vr );
static void init_s( lapack_int size, double *s );
static void init_sep( lapack_int size, double *sep );
static void init_work( lapack_int size, lapack_complex_double *work );
static void init_rwork( lapack_int size, double *rwork );
static int compare_ztrsna( double *s, double *s_i, double *sep, double *sep_i,
                           lapack_int m, lapack_int m_i, lapack_int info,
                           lapack_int info_i, char job, lapack_int mm );

int main(void)
{
    /* Local scalars */
    char job, job_i;
    char howmny, howmny_i;
    lapack_int n, n_i;
    lapack_int ldt, ldt_i;
    lapack_int ldt_r;
    lapack_int ldvl, ldvl_i;
    lapack_int ldvl_r;
    lapack_int ldvr, ldvr_i;
    lapack_int ldvr_r;
    lapack_int mm, mm_i;
    lapack_int m, m_i;
    lapack_int ldwork, ldwork_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_int *select = NULL, *select_i = NULL;
    lapack_complex_double *t = NULL, *t_i = NULL;
    lapack_complex_double *vl = NULL, *vl_i = NULL;
    lapack_complex_double *vr = NULL, *vr_i = NULL;
    double *s = NULL, *s_i = NULL;
    double *sep = NULL, *sep_i = NULL;
    lapack_complex_double *work = NULL, *work_i = NULL;
    double *rwork = NULL, *rwork_i = NULL;
    double *s_save = NULL;
    double *sep_save = NULL;
    lapack_complex_double *t_r = NULL;
    lapack_complex_double *vl_r = NULL;
    lapack_complex_double *vr_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ztrsna( &job, &howmny, &n, &ldt, &ldvl, &ldvr, &mm, &ldwork );
    ldt_r = n+2;
    ldvl_r = mm+2;
    ldvr_r = mm+2;
    job_i = job;
    howmny_i = howmny;
    n_i = n;
    ldt_i = ldt;
    ldvl_i = ldvl;
    ldvr_i = ldvr;
    mm_i = mm;
    ldwork_i = ldwork;

    /* Allocate memory for the LAPACK routine arrays */
    select = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t = (lapack_complex_double *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_double) );
    vl = (lapack_complex_double *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_double) );
    vr = (lapack_complex_double *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_double) );
    s = (double *)LAPACKE_malloc( mm * sizeof(double) );
    sep = (double *)LAPACKE_malloc( mm * sizeof(double) );
    work = (lapack_complex_double *)
        LAPACKE_malloc( ldwork*(n+1) * sizeof(lapack_complex_double) );
    rwork = (double *)LAPACKE_malloc( n * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t_i = (lapack_complex_double *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_double) );
    vl_i = (lapack_complex_double *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_double) );
    vr_i = (lapack_complex_double *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_double) );
    s_i = (double *)LAPACKE_malloc( mm * sizeof(double) );
    sep_i = (double *)LAPACKE_malloc( mm * sizeof(double) );
    work_i = (lapack_complex_double *)
        LAPACKE_malloc( ldwork*(n+1) * sizeof(lapack_complex_double) );
    rwork_i = (double *)LAPACKE_malloc( n * sizeof(double) );

    /* Allocate memory for the backup arrays */
    s_save = (double *)LAPACKE_malloc( mm * sizeof(double) );
    sep_save = (double *)LAPACKE_malloc( mm * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    t_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );
    vl_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(mm+2) * sizeof(lapack_complex_double) );
    vr_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(mm+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_select( n, select );
    init_t( ldt*n, t );
    init_vl( ldvl*mm, vl );
    init_vr( ldvr*mm, vr );
    init_s( mm, s );
    init_sep( mm, sep );
    init_work( ldwork*(n+1), work );
    init_rwork( n, rwork );

    /* Backup the ouptut arrays */
    for( i = 0; i < mm; i++ ) {
        s_save[i] = s[i];
    }
    for( i = 0; i < mm; i++ ) {
        sep_save[i] = sep[i];
    }

    /* Call the LAPACK routine */
    ztrsna_( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep,
             &mm, &m, work, &ldwork, rwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr[i];
    }
    for( i = 0; i < mm; i++ ) {
        s_i[i] = s_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        sep_i[i] = sep_save[i];
    }
    for( i = 0; i < ldwork*(n+1); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_ztrsna_work( LAPACK_COL_MAJOR, job_i, howmny_i, select_i,
                                  n_i, t_i, ldt_i, vl_i, ldvl_i, vr_i, ldvr_i,
                                  s_i, sep_i, mm_i, &m_i, work_i, ldwork_i,
                                  rwork_i );

    failed = compare_ztrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ztrsna\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ztrsna\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr[i];
    }
    for( i = 0; i < mm; i++ ) {
        s_i[i] = s_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        sep_i[i] = sep_save[i];
    }
    for( i = 0; i < ldwork*(n+1); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_ztrsna( LAPACK_COL_MAJOR, job_i, howmny_i, select_i, n_i,
                             t_i, ldt_i, vl_i, ldvl_i, vr_i, ldvr_i, s_i, sep_i,
                             mm_i, &m_i );

    failed = compare_ztrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ztrsna\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ztrsna\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr[i];
    }
    for( i = 0; i < mm; i++ ) {
        s_i[i] = s_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        sep_i[i] = sep_save[i];
    }
    for( i = 0; i < ldwork*(n+1); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_ztrsna_work( LAPACK_ROW_MAJOR, job_i, howmny_i, select_i,
                                  n_i, t_r, ldt_r, vl_r, ldvl_r, vr_r, ldvr_r,
                                  s_i, sep_i, mm_i, &m_i, work_i, ldwork_i,
                                  rwork_i );

    failed = compare_ztrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ztrsna\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ztrsna\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr[i];
    }
    for( i = 0; i < mm; i++ ) {
        s_i[i] = s_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        sep_i[i] = sep_save[i];
    }
    for( i = 0; i < ldwork*(n+1); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_ztrsna( LAPACK_ROW_MAJOR, job_i, howmny_i, select_i, n_i,
                             t_r, ldt_r, vl_r, ldvl_r, vr_r, ldvr_r, s_i, sep_i,
                             mm_i, &m_i );

    failed = compare_ztrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ztrsna\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ztrsna\n" );
    }

    /* Release memory */
    if( select != NULL ) {
        LAPACKE_free( select );
    }
    if( select_i != NULL ) {
        LAPACKE_free( select_i );
    }
    if( t != NULL ) {
        LAPACKE_free( t );
    }
    if( t_i != NULL ) {
        LAPACKE_free( t_i );
    }
    if( t_r != NULL ) {
        LAPACKE_free( t_r );
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
    if( vr != NULL ) {
        LAPACKE_free( vr );
    }
    if( vr_i != NULL ) {
        LAPACKE_free( vr_i );
    }
    if( vr_r != NULL ) {
        LAPACKE_free( vr_r );
    }
    if( s != NULL ) {
        LAPACKE_free( s );
    }
    if( s_i != NULL ) {
        LAPACKE_free( s_i );
    }
    if( s_save != NULL ) {
        LAPACKE_free( s_save );
    }
    if( sep != NULL ) {
        LAPACKE_free( sep );
    }
    if( sep_i != NULL ) {
        LAPACKE_free( sep_i );
    }
    if( sep_save != NULL ) {
        LAPACKE_free( sep_save );
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

/* Auxiliary function: ztrsna scalar parameters initialization */
static void init_scalars_ztrsna( char *job, char *howmny, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldvl,
                                 lapack_int *ldvr, lapack_int *mm,
                                 lapack_int *ldwork )
{
    *job = 'B';
    *howmny = 'A';
    *n = 4;
    *ldt = 8;
    *ldvl = 8;
    *ldvr = 8;
    *mm = 4;
    *ldwork = 8;

    return;
}

/* Auxiliary functions: ztrsna array parameters initialization */
static void init_select( lapack_int size, lapack_int *select ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        select[i] = 0;
    }
    select[0] = 0;
    select[1] = 0;
    select[2] = 0;
    select[3] = 0;
}
static void init_t( lapack_int size, lapack_complex_double *t ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        t[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    t[0] = lapack_make_complex_double( -6.00040000000000000e+000,
                                       -6.99990000000000020e+000 );
    t[8] = lapack_make_complex_double( 3.63700000000000020e-001,
                                       -3.65599999999999980e-001 );
    t[16] = lapack_make_complex_double( -1.88000000000000000e-001,
                                        4.78700000000000010e-001 );
    t[24] = lapack_make_complex_double( 8.78499999999999950e-001,
                                        -2.53900000000000010e-001 );
    t[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    t[9] = lapack_make_complex_double( -5.00000000000000000e+000,
                                       2.00599999999999980e+000 );
    t[17] = lapack_make_complex_double( -3.07000000000000020e-002,
                                        -7.21700000000000010e-001 );
    t[25] = lapack_make_complex_double( -2.29000000000000010e-001,
                                        1.31300000000000000e-001 );
    t[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    t[10] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    t[18] = lapack_make_complex_double( 7.99819999999999980e+000,
                                        -9.96399999999999950e-001 );
    t[26] = lapack_make_complex_double( 9.35699999999999980e-001,
                                        5.35900000000000040e-001 );
    t[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    t[11] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    t[19] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    t[27] = lapack_make_complex_double( 3.00230000000000000e+000,
                                        -3.99980000000000000e+000 );
}
static void init_vl( lapack_int size, lapack_complex_double *vl ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vl[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    vl[0] = lapack_make_complex_double( 1.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[8] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[16] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[24] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[1] = lapack_make_complex_double( 3.56694351594852160e-002,
                                        -4.43468951391364570e-002 );
    vl[9] = lapack_make_complex_double( 1.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vl[17] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[25] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[2] = lapack_make_complex_double( -2.20737607420687770e-003,
                                        3.13134125343339370e-002 );
    vl[10] = lapack_make_complex_double( -9.93319711341404500e-003,
                                         -5.32286446574668500e-002 );
    vl[18] = lapack_make_complex_double( 1.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[26] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vl[3] = lapack_make_complex_double( -7.82438806323640660e-002,
                                        -5.82707855845602950e-002 );
    vl[11] = lapack_make_complex_double( 3.18749533045133910e-002,
                                         -1.95590668724408750e-003 );
    vl[19] = lapack_make_complex_double( 1.84940888986473510e-001,
                                         3.91350226825490630e-003 );
    vl[27] = lapack_make_complex_double( 1.00000000000000000e+000,
                                         0.00000000000000000e+000 );
}
static void init_vr( lapack_int size, lapack_complex_double *vr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vr[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    vr[0] = lapack_make_complex_double( 1.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[8] = lapack_make_complex_double( -3.56694351594852160e-002,
                                        -4.43468951391364570e-002 );
    vr[16] = lapack_make_complex_double( -5.07460579179468000e-004,
                                         3.27715417727857880e-002 );
    vr[24] = lapack_make_complex_double( 7.92597025312869050e-002,
                                         -6.28502483030855260e-002 );
    vr[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[9] = lapack_make_complex_double( 1.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[17] = lapack_make_complex_double( 9.93319711341404500e-003,
                                         -5.32286446574668500e-002 );
    vr[25] = lapack_make_complex_double( -3.35036971875429250e-002,
                                         7.92711976468730460e-003 );
    vr[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[10] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[18] = lapack_make_complex_double( 1.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[26] = lapack_make_complex_double( -1.84940888986473510e-001,
                                         3.91350226825490630e-003 );
    vr[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vr[11] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[19] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vr[27] = lapack_make_complex_double( 1.00000000000000000e+000,
                                         0.00000000000000000e+000 );
}
static void init_s( lapack_int size, double *s ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        s[i] = 0;
    }
}
static void init_sep( lapack_int size, double *sep ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        sep[i] = 0;
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

/* Auxiliary function: C interface to ztrsna results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ztrsna( double *s, double *s_i, double *sep, double *sep_i,
                           lapack_int m, lapack_int m_i, lapack_int info,
                           lapack_int info_i, char job, lapack_int mm )
{
    lapack_int i;
    int failed = 0;
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        for( i = 0; i < mm; i++ ) {
            failed += compare_doubles(s[i],s_i[i]);
        }
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'v' ) ) {
        for( i = 0; i < mm; i++ ) {
            failed += compare_doubles(sep[i],sep_i[i]);
        }
    }
    failed += (m == m_i) ? 0 : 1;
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
