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
* dtrexc_1 is the test program for the C interface to LAPACK
* routine dtrexc
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

static void init_scalars_dtrexc( char *compq, lapack_int *n, lapack_int *ldt,
                                 lapack_int *ldq, lapack_int *ifst,
                                 lapack_int *ilst );
static void init_t( lapack_int size, double *t );
static void init_q( lapack_int size, double *q );
static void init_work( lapack_int size, double *work );
static int compare_dtrexc( double *t, double *t_i, double *q, double *q_i,
                           lapack_int ifst, lapack_int ifst_i, lapack_int ilst,
                           lapack_int ilst_i, lapack_int info,
                           lapack_int info_i, char compq, lapack_int ldq,
                           lapack_int ldt, lapack_int n );

int main(void)
{
    /* Local scalars */
    char compq, compq_i;
    lapack_int n, n_i;
    lapack_int ldt, ldt_i;
    lapack_int ldt_r;
    lapack_int ldq, ldq_i;
    lapack_int ldq_r;
    lapack_int ifst, ifst_i, ifst_save;
    lapack_int ilst, ilst_i, ilst_save;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *t = NULL, *t_i = NULL;
    double *q = NULL, *q_i = NULL;
    double *work = NULL, *work_i = NULL;
    double *t_save = NULL;
    double *q_save = NULL;
    double *t_r = NULL;
    double *q_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dtrexc( &compq, &n, &ldt, &ldq, &ifst, &ilst );
    ldt_r = n+2;
    ldq_r = n+2;
    compq_i = compq;
    n_i = n;
    ldt_i = ldt;
    ldq_i = ldq;
    ifst_i = ifst_save = ifst;
    ilst_i = ilst_save = ilst;

    /* Allocate memory for the LAPACK routine arrays */
    t = (double *)LAPACKE_malloc( ldt*n * sizeof(double) );
    q = (double *)LAPACKE_malloc( ldq*n * sizeof(double) );
    work = (double *)LAPACKE_malloc( n * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    t_i = (double *)LAPACKE_malloc( ldt*n * sizeof(double) );
    q_i = (double *)LAPACKE_malloc( ldq*n * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( n * sizeof(double) );

    /* Allocate memory for the backup arrays */
    t_save = (double *)LAPACKE_malloc( ldt*n * sizeof(double) );
    q_save = (double *)LAPACKE_malloc( ldq*n * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    t_r = (double *)LAPACKE_malloc( n*(n+2) * sizeof(double) );
    q_r = (double *)LAPACKE_malloc( n*(n+2) * sizeof(double) );

    /* Initialize input arrays */
    init_t( ldt*n, t );
    init_q( ldq*n, q );
    init_work( n, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldt*n; i++ ) {
        t_save[i] = t[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_save[i] = q[i];
    }

    /* Call the LAPACK routine */
    dtrexc_( &compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, work, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    ifst_i = ifst_save;
    ilst_i = ilst_save;
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dtrexc_work( LAPACK_COL_MAJOR, compq_i, n_i, t_i, ldt_i,
                                  q_i, ldq_i, &ifst_i, &ilst_i, work_i );

    failed = compare_dtrexc( t, t_i, q, q_i, ifst, ifst_i, ilst, ilst_i, info,
                             info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dtrexc\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dtrexc\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    ifst_i = ifst_save;
    ilst_i = ilst_save;
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dtrexc( LAPACK_COL_MAJOR, compq_i, n_i, t_i, ldt_i, q_i,
                             ldq_i, &ifst_i, &ilst_i );

    failed = compare_dtrexc( t, t_i, q, q_i, ifst, ifst_i, ilst, ilst_i, info,
                             info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dtrexc\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dtrexc\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    ifst_i = ifst_save;
    ilst_i = ilst_save;
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_dtrexc_work( LAPACK_ROW_MAJOR, compq_i, n_i, t_r, ldt_r,
                                  q_r, ldq_r, &ifst_i, &ilst_i, work_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, n, t_r, n+2, t_i, ldt );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_dtrexc( t, t_i, q, q_i, ifst, ifst_i, ilst, ilst_i, info,
                             info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dtrexc\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dtrexc\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    ifst_i = ifst_save;
    ilst_i = ilst_save;
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_dtrexc( LAPACK_ROW_MAJOR, compq_i, n_i, t_r, ldt_r, q_r,
                             ldq_r, &ifst_i, &ilst_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, n, t_r, n+2, t_i, ldt );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_dtrexc( t, t_i, q, q_i, ifst, ifst_i, ilst, ilst_i, info,
                             info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dtrexc\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dtrexc\n" );
    }

    /* Release memory */
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
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: dtrexc scalar parameters initialization */
static void init_scalars_dtrexc( char *compq, lapack_int *n, lapack_int *ldt,
                                 lapack_int *ldq, lapack_int *ifst,
                                 lapack_int *ilst )
{
    *compq = 'N';
    *n = 4;
    *ldt = 8;
    *ldq = 1;
    *ifst = 2;
    *ilst = 1;

    return;
}

/* Auxiliary functions: dtrexc array parameters initialization */
static void init_t( lapack_int size, double *t ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        t[i] = 0;
    }
    t[0] = 8.00000000000000040e-001;  /* t[0,0] */
    t[8] = -1.10000000000000000e-001;  /* t[0,1] */
    t[16] = 1.00000000000000000e-002;  /* t[0,2] */
    t[24] = 2.99999999999999990e-002;  /* t[0,3] */
    t[1] = 0.00000000000000000e+000;  /* t[1,0] */
    t[9] = -1.00000000000000010e-001;  /* t[1,1] */
    t[17] = 2.50000000000000000e-001;  /* t[1,2] */
    t[25] = 3.49999999999999980e-001;  /* t[1,3] */
    t[2] = 0.00000000000000000e+000;  /* t[2,0] */
    t[10] = -6.50000000000000020e-001;  /* t[2,1] */
    t[18] = -1.00000000000000010e-001;  /* t[2,2] */
    t[26] = 2.00000000000000010e-001;  /* t[2,3] */
    t[3] = 0.00000000000000000e+000;  /* t[3,0] */
    t[11] = 0.00000000000000000e+000;  /* t[3,1] */
    t[19] = 0.00000000000000000e+000;  /* t[3,2] */
    t[27] = -1.00000000000000010e-001;  /* t[3,3] */
}
static void init_q( lapack_int size, double *q ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        q[i] = 0;
    }
    q[0] = 0.00000000000000000e+000;  /* q[0,0] */
    q[1] = 0.00000000000000000e+000;  /* q[0,1] */
    q[2] = 0.00000000000000000e+000;  /* q[0,2] */
    q[3] = 0.00000000000000000e+000;  /* q[0,3] */
    q[1] = 0.00000000000000000e+000;  /* q[1,0] */
    q[2] = 0.00000000000000000e+000;  /* q[1,1] */
    q[3] = 0.00000000000000000e+000;  /* q[1,2] */
    q[2] = 0.00000000000000000e+000;  /* q[2,0] */
    q[3] = 0.00000000000000000e+000;  /* q[2,1] */
    q[3] = 0.00000000000000000e+000;  /* q[3,0] */
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to dtrexc results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dtrexc( double *t, double *t_i, double *q, double *q_i,
                           lapack_int ifst, lapack_int ifst_i, lapack_int ilst,
                           lapack_int ilst_i, lapack_int info,
                           lapack_int info_i, char compq, lapack_int ldq,
                           lapack_int ldt, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldt*n; i++ ) {
        failed += compare_doubles(t[i],t_i[i]);
    }
    if( LAPACKE_lsame( compq, 'v' ) ) {
        for( i = 0; i < ldq*n; i++ ) {
            failed += compare_doubles(q[i],q_i[i]);
        }
    }
    failed += (ifst == ifst_i) ? 0 : 1;
    failed += (ilst == ilst_i) ? 0 : 1;
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
