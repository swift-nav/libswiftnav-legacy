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
* ctrsna_1 is the test program for the C interface to LAPACK
* routine ctrsna
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

static void init_scalars_ctrsna( char *job, char *howmny, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldvl,
                                 lapack_int *ldvr, lapack_int *mm,
                                 lapack_int *ldwork );
static void init_select( lapack_int size, lapack_int *select );
static void init_t( lapack_int size, lapack_complex_float *t );
static void init_vl( lapack_int size, lapack_complex_float *vl );
static void init_vr( lapack_int size, lapack_complex_float *vr );
static void init_s( lapack_int size, float *s );
static void init_sep( lapack_int size, float *sep );
static void init_work( lapack_int size, lapack_complex_float *work );
static void init_rwork( lapack_int size, float *rwork );
static int compare_ctrsna( float *s, float *s_i, float *sep, float *sep_i,
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
    lapack_complex_float *t = NULL, *t_i = NULL;
    lapack_complex_float *vl = NULL, *vl_i = NULL;
    lapack_complex_float *vr = NULL, *vr_i = NULL;
    float *s = NULL, *s_i = NULL;
    float *sep = NULL, *sep_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    float *rwork = NULL, *rwork_i = NULL;
    float *s_save = NULL;
    float *sep_save = NULL;
    lapack_complex_float *t_r = NULL;
    lapack_complex_float *vl_r = NULL;
    lapack_complex_float *vr_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ctrsna( &job, &howmny, &n, &ldt, &ldvl, &ldvr, &mm, &ldwork );
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
    t = (lapack_complex_float *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_float) );
    vl = (lapack_complex_float *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_float) );
    vr = (lapack_complex_float *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_float) );
    s = (float *)LAPACKE_malloc( mm * sizeof(float) );
    sep = (float *)LAPACKE_malloc( mm * sizeof(float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( ldwork*(n+1) * sizeof(lapack_complex_float) );
    rwork = (float *)LAPACKE_malloc( n * sizeof(float) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t_i = (lapack_complex_float *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_float) );
    vl_i = (lapack_complex_float *)
        LAPACKE_malloc( ldvl*mm * sizeof(lapack_complex_float) );
    vr_i = (lapack_complex_float *)
        LAPACKE_malloc( ldvr*mm * sizeof(lapack_complex_float) );
    s_i = (float *)LAPACKE_malloc( mm * sizeof(float) );
    sep_i = (float *)LAPACKE_malloc( mm * sizeof(float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( ldwork*(n+1) * sizeof(lapack_complex_float) );
    rwork_i = (float *)LAPACKE_malloc( n * sizeof(float) );

    /* Allocate memory for the backup arrays */
    s_save = (float *)LAPACKE_malloc( mm * sizeof(float) );
    sep_save = (float *)LAPACKE_malloc( mm * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    t_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );
    vl_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(mm+2) * sizeof(lapack_complex_float) );
    vr_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(mm+2) * sizeof(lapack_complex_float) );

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
    ctrsna_( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep,
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
    info_i = LAPACKE_ctrsna_work( LAPACK_COL_MAJOR, job_i, howmny_i, select_i,
                                  n_i, t_i, ldt_i, vl_i, ldvl_i, vr_i, ldvr_i,
                                  s_i, sep_i, mm_i, &m_i, work_i, ldwork_i,
                                  rwork_i );

    failed = compare_ctrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ctrsna\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ctrsna\n" );
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
    info_i = LAPACKE_ctrsna( LAPACK_COL_MAJOR, job_i, howmny_i, select_i, n_i,
                             t_i, ldt_i, vl_i, ldvl_i, vr_i, ldvr_i, s_i, sep_i,
                             mm_i, &m_i );

    failed = compare_ctrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ctrsna\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ctrsna\n" );
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

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_ctrsna_work( LAPACK_ROW_MAJOR, job_i, howmny_i, select_i,
                                  n_i, t_r, ldt_r, vl_r, ldvl_r, vr_r, ldvr_r,
                                  s_i, sep_i, mm_i, &m_i, work_i, ldwork_i,
                                  rwork_i );

    failed = compare_ctrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ctrsna\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ctrsna\n" );
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
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_ctrsna( LAPACK_ROW_MAJOR, job_i, howmny_i, select_i, n_i,
                             t_r, ldt_r, vl_r, ldvl_r, vr_r, ldvr_r, s_i, sep_i,
                             mm_i, &m_i );

    failed = compare_ctrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ctrsna\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ctrsna\n" );
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

/* Auxiliary function: ctrsna scalar parameters initialization */
static void init_scalars_ctrsna( char *job, char *howmny, lapack_int *n,
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

/* Auxiliary functions: ctrsna array parameters initialization */
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
static void init_t( lapack_int size, lapack_complex_float *t ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        t[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    t[0] = lapack_make_complex_float( -6.000400066e+000, -6.999899864e+000 );
    t[8] = lapack_make_complex_float( 3.637000024e-001, -3.655999899e-001 );
    t[16] = lapack_make_complex_float( -1.879999936e-001, 4.787000120e-001 );
    t[24] = lapack_make_complex_float( 8.784999847e-001, -2.538999915e-001 );
    t[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[9] = lapack_make_complex_float( -5.000000000e+000, 2.006000042e+000 );
    t[17] = lapack_make_complex_float( -3.070000000e-002, -7.217000127e-001 );
    t[25] = lapack_make_complex_float( -2.290000021e-001, 1.313000023e-001 );
    t[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[10] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[18] = lapack_make_complex_float( 7.998199940e+000, -9.963999987e-001 );
    t[26] = lapack_make_complex_float( 9.356999993e-001, 5.358999968e-001 );
    t[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[19] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[27] = lapack_make_complex_float( 3.002300024e+000, -3.999799967e+000 );
}
static void init_vl( lapack_int size, lapack_complex_float *vl ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vl[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    vl[0] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
    vl[8] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[16] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[24] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[1] = lapack_make_complex_float( 3.566943482e-002, -4.434689879e-002 );
    vl[9] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
    vl[17] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[25] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[2] = lapack_make_complex_float( -2.207376994e-003, 3.131341189e-002 );
    vl[10] = lapack_make_complex_float( -9.933197871e-003, -5.322864652e-002 );
    vl[18] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
    vl[26] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vl[3] = lapack_make_complex_float( -7.824388146e-002, -5.827078223e-002 );
    vl[11] = lapack_make_complex_float( 3.187495098e-002, -1.955906162e-003 );
    vl[19] = lapack_make_complex_float( 1.849408895e-001, 3.913496155e-003 );
    vl[27] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
}
static void init_vr( lapack_int size, lapack_complex_float *vr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vr[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    vr[0] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
    vr[8] = lapack_make_complex_float( -3.566943482e-002, -4.434689879e-002 );
    vr[16] = lapack_make_complex_float( -5.074592773e-004, 3.277154267e-002 );
    vr[24] = lapack_make_complex_float( 7.925970852e-002, -6.285025179e-002 );
    vr[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[9] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
    vr[17] = lapack_make_complex_float( 9.933197871e-003, -5.322864652e-002 );
    vr[25] = lapack_make_complex_float( -3.350369632e-002, 7.927120663e-003 );
    vr[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[10] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[18] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
    vr[26] = lapack_make_complex_float( -1.849408895e-001, 3.913496155e-003 );
    vr[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[19] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    vr[27] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
}
static void init_s( lapack_int size, float *s ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        s[i] = 0;
    }
}
static void init_sep( lapack_int size, float *sep ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        sep[i] = 0;
    }
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

/* Auxiliary function: C interface to ctrsna results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ctrsna( float *s, float *s_i, float *sep, float *sep_i,
                           lapack_int m, lapack_int m_i, lapack_int info,
                           lapack_int info_i, char job, lapack_int mm )
{
    lapack_int i;
    int failed = 0;
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        for( i = 0; i < mm; i++ ) {
            failed += compare_floats(s[i],s_i[i]);
        }
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'v' ) ) {
        for( i = 0; i < mm; i++ ) {
            failed += compare_floats(sep[i],sep_i[i]);
        }
    }
    failed += (m == m_i) ? 0 : 1;
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
