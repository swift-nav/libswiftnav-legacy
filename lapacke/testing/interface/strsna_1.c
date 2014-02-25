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
* strsna_1 is the test program for the C interface to LAPACK
* routine strsna
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

static void init_scalars_strsna( char *job, char *howmny, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldvl,
                                 lapack_int *ldvr, lapack_int *mm,
                                 lapack_int *ldwork );
static void init_select( lapack_int size, lapack_int *select );
static void init_t( lapack_int size, float *t );
static void init_vl( lapack_int size, float *vl );
static void init_vr( lapack_int size, float *vr );
static void init_s( lapack_int size, float *s );
static void init_sep( lapack_int size, float *sep );
static void init_work( lapack_int size, float *work );
static void init_iwork( lapack_int size, lapack_int *iwork );
static int compare_strsna( float *s, float *s_i, float *sep, float *sep_i,
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
    float *t = NULL, *t_i = NULL;
    float *vl = NULL, *vl_i = NULL;
    float *vr = NULL, *vr_i = NULL;
    float *s = NULL, *s_i = NULL;
    float *sep = NULL, *sep_i = NULL;
    float *work = NULL, *work_i = NULL;
    lapack_int *iwork = NULL, *iwork_i = NULL;
    float *s_save = NULL;
    float *sep_save = NULL;
    float *t_r = NULL;
    float *vl_r = NULL;
    float *vr_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_strsna( &job, &howmny, &n, &ldt, &ldvl, &ldvr, &mm, &ldwork );
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
    t = (float *)LAPACKE_malloc( ldt*n * sizeof(float) );
    vl = (float *)LAPACKE_malloc( ldvl*mm * sizeof(float) );
    vr = (float *)LAPACKE_malloc( ldvr*mm * sizeof(float) );
    s = (float *)LAPACKE_malloc( mm * sizeof(float) );
    sep = (float *)LAPACKE_malloc( mm * sizeof(float) );
    work = (float *)LAPACKE_malloc( ldwork*(n+6) * sizeof(float) );
    iwork = (lapack_int *)LAPACKE_malloc( ((2*(n-1))) * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t_i = (float *)LAPACKE_malloc( ldt*n * sizeof(float) );
    vl_i = (float *)LAPACKE_malloc( ldvl*mm * sizeof(float) );
    vr_i = (float *)LAPACKE_malloc( ldvr*mm * sizeof(float) );
    s_i = (float *)LAPACKE_malloc( mm * sizeof(float) );
    sep_i = (float *)LAPACKE_malloc( mm * sizeof(float) );
    work_i = (float *)LAPACKE_malloc( ldwork*(n+6) * sizeof(float) );
    iwork_i = (lapack_int *)LAPACKE_malloc( ((2*(n-1))) * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    s_save = (float *)LAPACKE_malloc( mm * sizeof(float) );
    sep_save = (float *)LAPACKE_malloc( mm * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    t_r = (float *)LAPACKE_malloc( n*(n+2) * sizeof(float) );
    vl_r = (float *)LAPACKE_malloc( n*(mm+2) * sizeof(float) );
    vr_r = (float *)LAPACKE_malloc( n*(mm+2) * sizeof(float) );

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
    strsna_( &job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep,
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
    info_i = LAPACKE_strsna_work( LAPACK_COL_MAJOR, job_i, howmny_i, select_i,
                                  n_i, t_i, ldt_i, vl_i, ldvl_i, vr_i, ldvr_i,
                                  s_i, sep_i, mm_i, &m_i, work_i, ldwork_i,
                                  iwork_i );

    failed = compare_strsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to strsna\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to strsna\n" );
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
    info_i = LAPACKE_strsna( LAPACK_COL_MAJOR, job_i, howmny_i, select_i, n_i,
                             t_i, ldt_i, vl_i, ldvl_i, vr_i, ldvr_i, s_i, sep_i,
                             mm_i, &m_i );

    failed = compare_strsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to strsna\n" );
    } else {
        printf( "FAILED: column-major high-level interface to strsna\n" );
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

    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_strsna_work( LAPACK_ROW_MAJOR, job_i, howmny_i, select_i,
                                  n_i, t_r, ldt_r, vl_r, ldvl_r, vr_r, ldvr_r,
                                  s_i, sep_i, mm_i, &m_i, work_i, ldwork_i,
                                  iwork_i );

    failed = compare_strsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to strsna\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to strsna\n" );
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
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'e' ) ) {
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_strsna( LAPACK_ROW_MAJOR, job_i, howmny_i, select_i, n_i,
                             t_r, ldt_r, vl_r, ldvl_r, vr_r, ldvr_r, s_i, sep_i,
                             mm_i, &m_i );

    failed = compare_strsna( s, s_i, sep, sep_i, m, m_i, info, info_i, job,
                             mm );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to strsna\n" );
    } else {
        printf( "FAILED: row-major high-level interface to strsna\n" );
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

/* Auxiliary function: strsna scalar parameters initialization */
static void init_scalars_strsna( char *job, char *howmny, lapack_int *n,
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

/* Auxiliary functions: strsna array parameters initialization */
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
static void init_t( lapack_int size, float *t ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        t[i] = 0;
    }
    t[0] = 7.994999886e-001;  /* t[0,0] */
    t[8] = -1.143999994e-001;  /* t[0,1] */
    t[16] = 6.000000052e-003;  /* t[0,2] */
    t[24] = 3.359999880e-002;  /* t[0,3] */
    t[1] = 0.000000000e+000;  /* t[1,0] */
    t[9] = -9.939999878e-002;  /* t[1,1] */
    t[17] = 2.477999926e-001;  /* t[1,2] */
    t[25] = 3.474000096e-001;  /* t[1,3] */
    t[2] = 0.000000000e+000;  /* t[2,0] */
    t[10] = -6.482999921e-001;  /* t[2,1] */
    t[18] = -9.939999878e-002;  /* t[2,2] */
    t[26] = 2.026000023e-001;  /* t[2,3] */
    t[3] = 0.000000000e+000;  /* t[3,0] */
    t[11] = 0.000000000e+000;  /* t[3,1] */
    t[19] = 0.000000000e+000;  /* t[3,2] */
    t[27] = -1.006999984e-001;  /* t[3,3] */
}
static void init_vl( lapack_int size, float *vl ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vl[i] = 0;
    }
    vl[0] = 1.000000000e+000;  /* vl[0,0] */
    vl[8] = 0.000000000e+000;  /* vl[0,1] */
    vl[16] = 0.000000000e+000;  /* vl[0,2] */
    vl[24] = 0.000000000e+000;  /* vl[0,3] */
    vl[1] = -1.101757810e-001;  /* vl[1,0] */
    vl[9] = 8.492970467e-001;  /* vl[1,1] */
    vl[17] = 0.000000000e+000;  /* vl[1,2] */
    vl[25] = 0.000000000e+000;  /* vl[1,3] */
    vl[2] = -2.369735949e-002;  /* vl[2,0] */
    vl[10] = 0.000000000e+000;  /* vl[2,1] */
    vl[18] = 5.250760913e-001;  /* vl[2,2] */
    vl[26] = 0.000000000e+000;  /* vl[2,3] */
    vl[3] = -1.052671764e-002;  /* vl[3,0] */
    vl[11] = -2.630231678e-001;  /* vl[3,1] */
    vl[19] = 7.369768023e-001;  /* vl[3,2] */
    vl[27] = 1.000000000e+000;  /* vl[3,3] */
}
static void init_vr( lapack_int size, float *vr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vr[i] = 0;
    }
    vr[0] = 1.000000000e+000;  /* vr[0,0] */
    vr[8] = 6.811594218e-002;  /* vr[0,1] */
    vr[16] = 2.369735949e-002;  /* vr[0,2] */
    vr[24] = 8.112866431e-003;  /* vr[0,3] */
    vr[1] = 0.000000000e+000;  /* vr[1,0] */
    vr[9] = 6.182478666e-001;  /* vr[1,1] */
    vr[17] = 0.000000000e+000;  /* vr[1,2] */
    vr[25] = 2.206494659e-001;  /* vr[1,3] */
    vr[2] = 0.000000000e+000;  /* vr[2,0] */
    vr[10] = 0.000000000e+000;  /* vr[2,1] */
    vr[18] = 1.000000000e+000;  /* vr[2,2] */
    vr[26] = -1.000000000e+000;  /* vr[2,3] */
    vr[3] = 0.000000000e+000;  /* vr[3,0] */
    vr[11] = 0.000000000e+000;  /* vr[3,1] */
    vr[19] = 0.000000000e+000;  /* vr[3,2] */
    vr[27] = 7.124730945e-001;  /* vr[3,3] */
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
static void init_work( lapack_int size, float *work ) {
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

/* Auxiliary function: C interface to strsna results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_strsna( float *s, float *s_i, float *sep, float *sep_i,
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
