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
* ctrsen_1 is the test program for the C interface to LAPACK
* routine ctrsen
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

static void init_scalars_ctrsen( char *job, char *compq, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldq,
                                 lapack_int *lwork );
static void init_select( lapack_int size, lapack_int *select );
static void init_t( lapack_int size, lapack_complex_float *t );
static void init_q( lapack_int size, lapack_complex_float *q );
static void init_w( lapack_int size, lapack_complex_float *w );
static void init_work( lapack_int size, lapack_complex_float *work );
static int compare_ctrsen( lapack_complex_float *t, lapack_complex_float *t_i,
                           lapack_complex_float *q, lapack_complex_float *q_i,
                           lapack_complex_float *w, lapack_complex_float *w_i,
                           lapack_int m, lapack_int m_i, float s, float s_i,
                           float sep, float sep_i, lapack_int info,
                           lapack_int info_i, char compq, lapack_int ldq,
                           lapack_int ldt, lapack_int n );

int main(void)
{
    /* Local scalars */
    char job, job_i;
    char compq, compq_i;
    lapack_int n, n_i;
    lapack_int ldt, ldt_i;
    lapack_int ldt_r;
    lapack_int ldq, ldq_i;
    lapack_int ldq_r;
    lapack_int m, m_i;
    float s, s_i;
    float sep, sep_i;
    lapack_int lwork, lwork_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_int *select = NULL, *select_i = NULL;
    lapack_complex_float *t = NULL, *t_i = NULL;
    lapack_complex_float *q = NULL, *q_i = NULL;
    lapack_complex_float *w = NULL, *w_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    lapack_complex_float *t_save = NULL;
    lapack_complex_float *q_save = NULL;
    lapack_complex_float *w_save = NULL;
    lapack_complex_float *t_r = NULL;
    lapack_complex_float *q_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ctrsen( &job, &compq, &n, &ldt, &ldq, &lwork );
    ldt_r = n+2;
    ldq_r = n+2;
    job_i = job;
    compq_i = compq;
    n_i = n;
    ldt_i = ldt;
    ldq_i = ldq;
    lwork_i = lwork;

    /* Allocate memory for the LAPACK routine arrays */
    select = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t = (lapack_complex_float *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_float) );
    q = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );
    w = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t_i = (lapack_complex_float *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_float) );
    q_i = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );
    w_i = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    t_save = (lapack_complex_float *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_float) );
    q_save = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );
    w_save = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    t_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );
    q_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_select( n, select );
    init_t( ldt*n, t );
    init_q( ldq*n, q );
    init_w( n, w );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldt*n; i++ ) {
        t_save[i] = t[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_save[i] = q[i];
    }
    for( i = 0; i < n; i++ ) {
        w_save[i] = w[i];
    }

    /* Call the LAPACK routine */
    ctrsen_( &job, &compq, select, &n, t, &ldt, q, &ldq, w, &m, &s, &sep, work,
             &lwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_ctrsen_work( LAPACK_COL_MAJOR, job_i, compq_i, select_i,
                                  n_i, t_i, ldt_i, q_i, ldq_i, w_i, &m_i, &s_i,
                                  &sep_i, work_i, lwork_i );

    failed = compare_ctrsen( t, t_i, q, q_i, w, w_i, m, m_i, s, s_i, sep, sep_i,
                             info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ctrsen\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ctrsen\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_ctrsen( LAPACK_COL_MAJOR, job_i, compq_i, select_i, n_i,
                             t_i, ldt_i, q_i, ldq_i, w_i, &m_i, &s_i, &sep_i );

    failed = compare_ctrsen( t, t_i, q, q_i, w, w_i, m, m_i, s, s_i, sep, sep_i,
                             info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ctrsen\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ctrsen\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_ctrsen_work( LAPACK_ROW_MAJOR, job_i, compq_i, select_i,
                                  n_i, t_r, ldt_r, q_r, ldq_r, w_i, &m_i, &s_i,
                                  &sep_i, work_i, lwork_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, t_r, n+2, t_i, ldt );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_ctrsen( t, t_i, q, q_i, w, w_i, m, m_i, s, s_i, sep, sep_i,
                             info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ctrsen\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ctrsen\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select[i];
    }
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_ctrsen( LAPACK_ROW_MAJOR, job_i, compq_i, select_i, n_i,
                             t_r, ldt_r, q_r, ldq_r, w_i, &m_i, &s_i, &sep_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, t_r, n+2, t_i, ldt );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_ctrsen( t, t_i, q, q_i, w, w_i, m, m_i, s, s_i, sep, sep_i,
                             info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ctrsen\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ctrsen\n" );
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
    if( t_save != NULL ) {
        LAPACKE_free( t_save );
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
    if( w != NULL ) {
        LAPACKE_free( w );
    }
    if( w_i != NULL ) {
        LAPACKE_free( w_i );
    }
    if( w_save != NULL ) {
        LAPACKE_free( w_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: ctrsen scalar parameters initialization */
static void init_scalars_ctrsen( char *job, char *compq, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldq,
                                 lapack_int *lwork )
{
    *job = 'B';
    *compq = 'V';
    *n = 4;
    *ldt = 8;
    *ldq = 8;
    *lwork = 32;

    return;
}

/* Auxiliary functions: ctrsen array parameters initialization */
static void init_select( lapack_int size, lapack_int *select ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        select[i] = 0;
    }
    select[0] = -1;
    select[1] = 0;
    select[2] = 0;
    select[3] = -1;
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
static void init_q( lapack_int size, lapack_complex_float *q ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        q[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    q[0] = lapack_make_complex_float( -8.346999884e-001, -1.363999993e-001 );
    q[8] = lapack_make_complex_float( -6.279999763e-002, 3.806000054e-001 );
    q[16] = lapack_make_complex_float( 2.764999866e-001, -8.460000157e-002 );
    q[24] = lapack_make_complex_float( 6.329999864e-002, -2.198999971e-001 );
    q[1] = lapack_make_complex_float( 6.639999896e-002, -2.967999876e-001 );
    q[9] = lapack_make_complex_float( 2.364999950e-001, 5.239999890e-001 );
    q[17] = lapack_make_complex_float( -5.877000093e-001, -4.208000004e-001 );
    q[25] = lapack_make_complex_float( 8.349999785e-002, 2.182999998e-001 );
    q[2] = lapack_make_complex_float( -3.620000184e-002, -3.215000033e-001 );
    q[10] = lapack_make_complex_float( 3.143000007e-001, -5.472999811e-001 );
    q[18] = lapack_make_complex_float( 5.759999901e-002, -5.735999942e-001 );
    q[26] = lapack_make_complex_float( 5.700000096e-003, -4.058000147e-001 );
    q[3] = lapack_make_complex_float( 8.600000292e-003, 2.958000004e-001 );
    q[11] = lapack_make_complex_float( -3.416000009e-001, -7.569999993e-002 );
    q[19] = lapack_make_complex_float( -1.899999976e-001, -1.599999964e-001 );
    q[27] = lapack_make_complex_float( 8.327000141e-001, -1.868000031e-001 );
}
static void init_w( lapack_int size, lapack_complex_float *w ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        w[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}

/* Auxiliary function: C interface to ctrsen results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ctrsen( lapack_complex_float *t, lapack_complex_float *t_i,
                           lapack_complex_float *q, lapack_complex_float *q_i,
                           lapack_complex_float *w, lapack_complex_float *w_i,
                           lapack_int m, lapack_int m_i, float s, float s_i,
                           float sep, float sep_i, lapack_int info,
                           lapack_int info_i, char compq, lapack_int ldq,
                           lapack_int ldt, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldt*n; i++ ) {
        failed += compare_complex_floats(t[i],t_i[i]);
    }
    if( LAPACKE_lsame( compq, 'v' ) ) {
        for( i = 0; i < ldq*n; i++ ) {
            failed += compare_complex_floats(q[i],q_i[i]);
        }
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_complex_floats(w[i],w_i[i]);
    }
    failed += (m == m_i) ? 0 : 1;
    failed += compare_floats(s,s_i);
    failed += compare_floats(sep,sep_i);
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
