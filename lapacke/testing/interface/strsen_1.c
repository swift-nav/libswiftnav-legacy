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
* strsen_1 is the test program for the C interface to LAPACK
* routine strsen
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

static void init_scalars_strsen( char *job, char *compq, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldq,
                                 lapack_int *lwork, lapack_int *liwork );
static void init_select( lapack_int size, lapack_int *select );
static void init_t( lapack_int size, float *t );
static void init_q( lapack_int size, float *q );
static void init_wr( lapack_int size, float *wr );
static void init_wi( lapack_int size, float *wi );
static void init_work( lapack_int size, float *work );
static void init_iwork( lapack_int size, lapack_int *iwork );
static int compare_strsen( float *t, float *t_i, float *q, float *q_i,
                           float *wr, float *wr_i, float *wi, float *wi_i,
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
    lapack_int liwork, liwork_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_int *select = NULL, *select_i = NULL;
    float *t = NULL, *t_i = NULL;
    float *q = NULL, *q_i = NULL;
    float *wr = NULL, *wr_i = NULL;
    float *wi = NULL, *wi_i = NULL;
    float *work = NULL, *work_i = NULL;
    lapack_int *iwork = NULL, *iwork_i = NULL;
    float *t_save = NULL;
    float *q_save = NULL;
    float *wr_save = NULL;
    float *wi_save = NULL;
    float *t_r = NULL;
    float *q_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_strsen( &job, &compq, &n, &ldt, &ldq, &lwork, &liwork );
    ldt_r = n+2;
    ldq_r = n+2;
    job_i = job;
    compq_i = compq;
    n_i = n;
    ldt_i = ldt;
    ldq_i = ldq;
    lwork_i = lwork;
    liwork_i = liwork;

    /* Allocate memory for the LAPACK routine arrays */
    select = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t = (float *)LAPACKE_malloc( ldt*n * sizeof(float) );
    q = (float *)LAPACKE_malloc( ldq*n * sizeof(float) );
    wr = (float *)LAPACKE_malloc( n * sizeof(float) );
    wi = (float *)LAPACKE_malloc( n * sizeof(float) );
    work = (float *)LAPACKE_malloc( lwork * sizeof(float) );
    iwork = (lapack_int *)LAPACKE_malloc( liwork * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    t_i = (float *)LAPACKE_malloc( ldt*n * sizeof(float) );
    q_i = (float *)LAPACKE_malloc( ldq*n * sizeof(float) );
    wr_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    wi_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    work_i = (float *)LAPACKE_malloc( lwork * sizeof(float) );
    iwork_i = (lapack_int *)LAPACKE_malloc( liwork * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    t_save = (float *)LAPACKE_malloc( ldt*n * sizeof(float) );
    q_save = (float *)LAPACKE_malloc( ldq*n * sizeof(float) );
    wr_save = (float *)LAPACKE_malloc( n * sizeof(float) );
    wi_save = (float *)LAPACKE_malloc( n * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    t_r = (float *)LAPACKE_malloc( n*(n+2) * sizeof(float) );
    q_r = (float *)LAPACKE_malloc( n*(n+2) * sizeof(float) );

    /* Initialize input arrays */
    init_select( n, select );
    init_t( ldt*n, t );
    init_q( ldq*n, q );
    init_wr( n, wr );
    init_wi( n, wi );
    init_work( lwork, work );
    init_iwork( liwork, iwork );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldt*n; i++ ) {
        t_save[i] = t[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_save[i] = q[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_save[i] = wr[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_save[i] = wi[i];
    }

    /* Call the LAPACK routine */
    strsen_( &job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, &m, &s, &sep,
             work, &lwork, iwork, &liwork, &info );

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
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < liwork; i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_strsen_work( LAPACK_COL_MAJOR, job_i, compq_i, select_i,
                                  n_i, t_i, ldt_i, q_i, ldq_i, wr_i, wi_i, &m_i,
                                  &s_i, &sep_i, work_i, lwork_i, iwork_i,
                                  liwork_i );

    failed = compare_strsen( t, t_i, q, q_i, wr, wr_i, wi, wi_i, m, m_i, s, s_i,
                             sep, sep_i, info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to strsen\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to strsen\n" );
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
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < liwork; i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_strsen( LAPACK_COL_MAJOR, job_i, compq_i, select_i, n_i,
                             t_i, ldt_i, q_i, ldq_i, wr_i, wi_i, &m_i, &s_i,
                             &sep_i );

    failed = compare_strsen( t, t_i, q, q_i, wr, wr_i, wi, wi_i, m, m_i, s, s_i,
                             sep, sep_i, info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to strsen\n" );
    } else {
        printf( "FAILED: column-major high-level interface to strsen\n" );
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
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < liwork; i++ ) {
        iwork_i[i] = iwork[i];
    }

    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_strsen_work( LAPACK_ROW_MAJOR, job_i, compq_i, select_i,
                                  n_i, t_r, ldt_r, q_r, ldq_r, wr_i, wi_i, &m_i,
                                  &s_i, &sep_i, work_i, lwork_i, iwork_i,
                                  liwork_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, n, t_r, n+2, t_i, ldt );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_strsen( t, t_i, q, q_i, wr, wr_i, wi, wi_i, m, m_i, s, s_i,
                             sep, sep_i, info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to strsen\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to strsen\n" );
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
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < liwork; i++ ) {
        iwork_i[i] = iwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_strsen( LAPACK_ROW_MAJOR, job_i, compq_i, select_i, n_i,
                             t_r, ldt_r, q_r, ldq_r, wr_i, wi_i, &m_i, &s_i,
                             &sep_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, n, t_r, n+2, t_i, ldt );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_strsen( t, t_i, q, q_i, wr, wr_i, wi, wi_i, m, m_i, s, s_i,
                             sep, sep_i, info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to strsen\n" );
    } else {
        printf( "FAILED: row-major high-level interface to strsen\n" );
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
    if( wr != NULL ) {
        LAPACKE_free( wr );
    }
    if( wr_i != NULL ) {
        LAPACKE_free( wr_i );
    }
    if( wr_save != NULL ) {
        LAPACKE_free( wr_save );
    }
    if( wi != NULL ) {
        LAPACKE_free( wi );
    }
    if( wi_i != NULL ) {
        LAPACKE_free( wi_i );
    }
    if( wi_save != NULL ) {
        LAPACKE_free( wi_save );
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

/* Auxiliary function: strsen scalar parameters initialization */
static void init_scalars_strsen( char *job, char *compq, lapack_int *n,
                                 lapack_int *ldt, lapack_int *ldq,
                                 lapack_int *lwork, lapack_int *liwork )
{
    *job = 'B';
    *compq = 'V';
    *n = 4;
    *ldt = 8;
    *ldq = 8;
    *lwork = 32;
    *liwork = 32;

    return;
}

/* Auxiliary functions: strsen array parameters initialization */
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
static void init_q( lapack_int size, float *q ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        q[i] = 0;
    }
    q[0] = 6.550999880e-001;  /* q[0,0] */
    q[8] = 1.036999971e-001;  /* q[0,1] */
    q[16] = 3.449999988e-001;  /* q[0,2] */
    q[24] = 6.640999913e-001;  /* q[0,3] */
    q[1] = 5.235999823e-001;  /* q[1,0] */
    q[9] = -5.806999803e-001;  /* q[1,1] */
    q[17] = -6.140999794e-001;  /* q[1,2] */
    q[25] = -1.067999974e-001;  /* q[1,3] */
    q[2] = -5.361999869e-001;  /* q[2,0] */
    q[10] = -3.073000014e-001;  /* q[2,1] */
    q[18] = -2.935000062e-001;  /* q[2,2] */
    q[26] = 7.293000221e-001;  /* q[2,3] */
    q[3] = 9.560000151e-002;  /* q[3,0] */
    q[11] = 7.466999888e-001;  /* q[3,1] */
    q[19] = -6.463000178e-001;  /* q[3,2] */
    q[27] = 1.248999983e-001;  /* q[3,3] */
}
static void init_wr( lapack_int size, float *wr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        wr[i] = 0;
    }
}
static void init_wi( lapack_int size, float *wi ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        wi[i] = 0;
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

/* Auxiliary function: C interface to strsen results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_strsen( float *t, float *t_i, float *q, float *q_i,
                           float *wr, float *wr_i, float *wi, float *wi_i,
                           lapack_int m, lapack_int m_i, float s, float s_i,
                           float sep, float sep_i, lapack_int info,
                           lapack_int info_i, char compq, lapack_int ldq,
                           lapack_int ldt, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldt*n; i++ ) {
        failed += compare_floats(t[i],t_i[i]);
    }
    if( LAPACKE_lsame( compq, 'v' ) ) {
        for( i = 0; i < ldq*n; i++ ) {
            failed += compare_floats(q[i],q_i[i]);
        }
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_floats(wr[i],wr_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_floats(wi[i],wi_i[i]);
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
