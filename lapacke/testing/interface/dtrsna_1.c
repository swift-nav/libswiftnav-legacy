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
* dtrsna_1 is the test program for the C interface to LAPACK
* routine dtrsna
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

static void init_scalars_dtrsna( char *job, char *howmny, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldvl,
                                 lapack_int *ldvr, lapack_int *mm,
                                 lapack_int *ldwork );
static void init_select( lapack_int size, lapack_int *select );
static void init_t( lapack_int size, double *t );
static void init_vl( lapack_int size, double *vl );
static void init_vr( lapack_int size, double *vr );
static void init_s( lapack_int size, double *s );
static void init_sep( lapack_int size, double *sep );
static void init_work( lapack_int size, double *work );
static void init_iwork( lapack_int size, lapack_int *iwork );
static int compare_dtrsna( double *s, double *s_i, double *sep, double *sep_i,
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
    double *t = NULL, *t_i = NULL;
    double *vl = NULL, *vl_i = NULL;
    double *vr = NULL, *vr_i = NULL;
    double *s = NULL, *s_i = NULL;
    double *sep = NULL, *sep_i = NULL;
    double *work = NULL, *work_i = NULL;
    lapack_int *iwork = NULL, *iwork_i = NULL;
    double *s_save = NULL;
    double *sep_save = NULL;
    double *t_r = NULL;
    double *vl_r = NULL;
    double *vr_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dtrsna( &job, &howmny, &n, &ldt, &ldvl, &ldvr, &mm, &ldwork );
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
    t = (double *)LAPACKE_malloc( ldt*n * sizeof(double) );
    vl = (double *)LAPACKE_malloc( ldvl*mm * sizeof(double) );
    vr = (double *)LAPACKE_malloc( ldvr*mm * sizeof(double) );
    s = (double *)LAPACKE_malloc( mm * sizeof(double) );
    sep = (double *)LAPACKE_malloc( mm * sizeof(double) );
    work = (double *)LAPACKE_malloc( ldwork*(n+6) * sizeof(double) );
    iwork = (lapack_int *)LAPACKE_malloc( ((2*(n-1))) * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t_i = (double *)LAPACKE_malloc( ldt*n * sizeof(double) );
    vl_i = (double *)LAPACKE_malloc( ldvl*mm * sizeof(double) );
    vr_i = (double *)LAPACKE_malloc( ldvr*mm * sizeof(double) );
    s_i = (double *)LAPACKE_malloc( mm * sizeof(double) );
    sep_i = (double *)LAPACKE_malloc( mm * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( ldwork*(n+6) * sizeof(double) );
    iwork_i = (lapack_int *)LAPACKE_malloc( ((2*(n-1))) * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    s_save = (double *)LAPACKE_malloc( mm * sizeof(double) );
    sep_save = (double *)LAPACKE_malloc( mm * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    t_r = (double *)LAPACKE_malloc( n*(n+2) * sizeof(double) );
    vl_r = (double *)LAPACKE_malloc( n*(mm+2) * sizeof(double) );
    vr_r = (double *)LAPACKE_malloc( n*(mm+2) * sizeof(double) );

    /* Initialize input arrays */
    init_select( n, select );
    init_t( ldt*n, t );
    init_vl( ldvl*mm, vl );
    init_vr( ldvr*mm, vr );
    init_s( mm, s );
    init_sep( mm, sep );
    init_work( ldwork*(n+6), work );
    init_iwork( (2*(n-1)), iwork );

    /* Backup the ouptut arrays */
    for( i = 0; i < mm; i++ ) {
        s_save[i] = s[i];
    }
    for( i = 0; i < mm; i++ ) {
        sep_save[i] = sep[i];
    }

    /* Call the LAPACK routine */
    dtrsna_( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep,
             &mm, &m, work, &ldwork, iwork, &info );

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
    for( i = 0; i < ldwork*(n+6); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < (2*(n-1)); i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_dtrsna_work( LAPACK_COL_MAJOR, job_i, howmny_i, select_i,
                                  n_i, t_i, ldt_i, vl_i, ldvl_i, vr_i, ldvr_i,
                                  s_i, sep_i, mm_i, &m_i, work_i, ldwork_i,
                                  iwork_i );

    failed = compare_dtrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dtrsna\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dtrsna\n" );
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
    for( i = 0; i < ldwork*(n+6); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < (2*(n-1)); i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_dtrsna( LAPACK_COL_MAJOR, job_i, howmny_i, select_i, n_i,
                             t_i, ldt_i, vl_i, ldvl_i, vr_i, ldvr_i, s_i, sep_i,
                             mm_i, &m_i );

    failed = compare_dtrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dtrsna\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dtrsna\n" );
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
    for( i = 0; i < ldwork*(n+6); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < (2*(n-1)); i++ ) {
        iwork_i[i] = iwork[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_dtrsna_work( LAPACK_ROW_MAJOR, job_i, howmny_i, select_i,
                                  n_i, t_r, ldt_r, vl_r, ldvl_r, vr_r, ldvr_r,
                                  s_i, sep_i, mm_i, &m_i, work_i, ldwork_i,
                                  iwork_i );

    failed = compare_dtrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dtrsna\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dtrsna\n" );
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
    for( i = 0; i < ldwork*(n+6); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < (2*(n-1)); i++ ) {
        iwork_i[i] = iwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_dtrsna( LAPACK_ROW_MAJOR, job_i, howmny_i, select_i, n_i,
                             t_r, ldt_r, vl_r, ldvl_r, vr_r, ldvr_r, s_i, sep_i,
                             mm_i, &m_i );

    failed = compare_dtrsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dtrsna\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dtrsna\n" );
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
    if( iwork != NULL ) {
        LAPACKE_free( iwork );
    }
    if( iwork_i != NULL ) {
        LAPACKE_free( iwork_i );
    }

    return 0;
}

/* Auxiliary function: dtrsna scalar parameters initialization */
static void init_scalars_dtrsna( char *job, char *howmny, lapack_int *n,
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

/* Auxiliary functions: dtrsna array parameters initialization */
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
static void init_t( lapack_int size, double *t ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        t[i] = 0;
    }
    t[0] = 7.99499999999999990e-001;  /* t[0,0] */
    t[8] = -1.14400000000000000e-001;  /* t[0,1] */
    t[16] = 6.00000000000000010e-003;  /* t[0,2] */
    t[24] = 3.35999999999999980e-002;  /* t[0,3] */
    t[1] = 0.00000000000000000e+000;  /* t[1,0] */
    t[9] = -9.94000000000000020e-002;  /* t[1,1] */
    t[17] = 2.47799999999999990e-001;  /* t[1,2] */
    t[25] = 3.47399999999999990e-001;  /* t[1,3] */
    t[2] = 0.00000000000000000e+000;  /* t[2,0] */
    t[10] = -6.48299999999999990e-001;  /* t[2,1] */
    t[18] = -9.94000000000000020e-002;  /* t[2,2] */
    t[26] = 2.02600000000000000e-001;  /* t[2,3] */
    t[3] = 0.00000000000000000e+000;  /* t[3,0] */
    t[11] = 0.00000000000000000e+000;  /* t[3,1] */
    t[19] = 0.00000000000000000e+000;  /* t[3,2] */
    t[27] = -1.00700000000000000e-001;  /* t[3,3] */
}
static void init_vl( lapack_int size, double *vl ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vl[i] = 0;
    }
    vl[0] = 1.00000000000000000e+000;  /* vl[0,0] */
    vl[8] = 0.00000000000000000e+000;  /* vl[0,1] */
    vl[16] = 0.00000000000000000e+000;  /* vl[0,2] */
    vl[24] = 0.00000000000000000e+000;  /* vl[0,3] */
    vl[1] = -1.10175772459959130e-001;  /* vl[1,0] */
    vl[9] = 8.49297113677700110e-001;  /* vl[1,1] */
    vl[17] = 0.00000000000000000e+000;  /* vl[1,2] */
    vl[25] = 0.00000000000000000e+000;  /* vl[1,3] */
    vl[2] = -2.36973594566446470e-002;  /* vl[2,0] */
    vl[10] = 0.00000000000000000e+000;  /* vl[2,1] */
    vl[18] = 5.25076145338855160e-001;  /* vl[2,2] */
    vl[26] = 0.00000000000000000e+000;  /* vl[2,3] */
    vl[3] = -1.05267144840102290e-002;  /* vl[3,0] */
    vl[11] = -2.63023202214821120e-001;  /* vl[3,1] */
    vl[19] = 7.36976797785178770e-001;  /* vl[3,2] */
    vl[27] = 1.00000000000000000e+000;  /* vl[3,3] */
}
static void init_vr( lapack_int size, double *vr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vr[i] = 0;
    }
    vr[0] = 1.00000000000000000e+000;  /* vr[0,0] */
    vr[8] = 6.81159384405489730e-002;  /* vr[0,1] */
    vr[16] = 2.36973594566446470e-002;  /* vr[0,2] */
    vr[24] = 8.11286711399037890e-003;  /* vr[0,3] */
    vr[1] = 0.00000000000000000e+000;  /* vr[1,0] */
    vr[9] = 6.18247886261057530e-001;  /* vr[1,1] */
    vr[17] = 0.00000000000000000e+000;  /* vr[1,2] */
    vr[25] = 2.20649468606918140e-001;  /* vr[1,3] */
    vr[2] = 0.00000000000000000e+000;  /* vr[2,0] */
    vr[10] = 0.00000000000000000e+000;  /* vr[2,1] */
    vr[18] = 1.00000000000000000e+000;  /* vr[2,2] */
    vr[26] = -1.00000000000000000e+000;  /* vr[2,3] */
    vr[3] = 0.00000000000000000e+000;  /* vr[3,0] */
    vr[11] = 0.00000000000000000e+000;  /* vr[3,1] */
    vr[19] = 0.00000000000000000e+000;  /* vr[3,2] */
    vr[27] = 7.12473102161229170e-001;  /* vr[3,3] */
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
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}
static void init_iwork( lapack_int size, lapack_int *iwork ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        iwork[i] = 0;
    }
}

/* Auxiliary function: C interface to dtrsna results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dtrsna( double *s, double *s_i, double *sep, double *sep_i,
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
