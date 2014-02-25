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
* dgeqpf_1 is the test program for the C interface to LAPACK
* routine dgeqpf
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

static void init_scalars_dgeqpf( lapack_int *m, lapack_int *n,
                                 lapack_int *lda );
static void init_a( lapack_int size, double *a );
static void init_jpvt( lapack_int size, lapack_int *jpvt );
static void init_tau( lapack_int size, double *tau );
static void init_work( lapack_int size, double *work );
static int compare_dgeqpf( double *a, double *a_i, lapack_int *jpvt,
                           lapack_int *jpvt_i, double *tau, double *tau_i,
                           lapack_int info, lapack_int info_i, lapack_int lda,
                           lapack_int m, lapack_int n );

int main(void)
{
    /* Local scalars */
    lapack_int m, m_i;
    lapack_int n, n_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *a = NULL, *a_i = NULL;
    lapack_int *jpvt = NULL, *jpvt_i = NULL;
    double *tau = NULL, *tau_i = NULL;
    double *work = NULL, *work_i = NULL;
    double *a_save = NULL;
    lapack_int *jpvt_save = NULL;
    double *tau_save = NULL;
    double *a_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dgeqpf( &m, &n, &lda );
    lda_r = n+2;
    m_i = m;
    n_i = n;
    lda_i = lda;

    /* Allocate memory for the LAPACK routine arrays */
    a = (double *)LAPACKE_malloc( lda*n * sizeof(double) );
    jpvt = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    tau = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );
    work = (double *)LAPACKE_malloc( 3*n * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    a_i = (double *)LAPACKE_malloc( lda*n * sizeof(double) );
    jpvt_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    tau_i = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( 3*n * sizeof(double) );

    /* Allocate memory for the backup arrays */
    a_save = (double *)LAPACKE_malloc( lda*n * sizeof(double) );
    jpvt_save = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    tau_save = (double *)LAPACKE_malloc( MIN(m,n) * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    a_r = (double *)LAPACKE_malloc( m*(n+2) * sizeof(double) );

    /* Initialize input arrays */
    init_a( lda*n, a );
    init_jpvt( n, jpvt );
    init_tau( (MIN(m,n)), tau );
    init_work( 3*n, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < lda*n; i++ ) {
        a_save[i] = a[i];
    }
    for( i = 0; i < n; i++ ) {
        jpvt_save[i] = jpvt[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_save[i] = tau[i];
    }

    /* Call the LAPACK routine */
    dgeqpf_( &m, &n, a, &lda, jpvt, tau, work, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < n; i++ ) {
        jpvt_i[i] = jpvt_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_i[i] = tau_save[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dgeqpf_work( LAPACK_COL_MAJOR, m_i, n_i, a_i, lda_i,
                                  jpvt_i, tau_i, work_i );

    failed = compare_dgeqpf( a, a_i, jpvt, jpvt_i, tau, tau_i, info, info_i,
                             lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dgeqpf\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dgeqpf\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < n; i++ ) {
        jpvt_i[i] = jpvt_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_i[i] = tau_save[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dgeqpf( LAPACK_COL_MAJOR, m_i, n_i, a_i, lda_i, jpvt_i,
                             tau_i );

    failed = compare_dgeqpf( a, a_i, jpvt, jpvt_i, tau, tau_i, info, info_i,
                             lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dgeqpf\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dgeqpf\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < n; i++ ) {
        jpvt_i[i] = jpvt_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_i[i] = tau_save[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_dgeqpf_work( LAPACK_ROW_MAJOR, m_i, n_i, a_r, lda_r,
                                  jpvt_i, tau_i, work_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, n, a_r, n+2, a_i, lda );

    failed = compare_dgeqpf( a, a_i, jpvt, jpvt_i, tau, tau_i, info, info_i,
                             lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dgeqpf\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dgeqpf\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < n; i++ ) {
        jpvt_i[i] = jpvt_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tau_i[i] = tau_save[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_dgeqpf( LAPACK_ROW_MAJOR, m_i, n_i, a_r, lda_r, jpvt_i,
                             tau_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, n, a_r, n+2, a_i, lda );

    failed = compare_dgeqpf( a, a_i, jpvt, jpvt_i, tau, tau_i, info, info_i,
                             lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dgeqpf\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dgeqpf\n" );
    }

    /* Release memory */
    if( a != NULL ) {
        LAPACKE_free( a );
    }
    if( a_i != NULL ) {
        LAPACKE_free( a_i );
    }
    if( a_r != NULL ) {
        LAPACKE_free( a_r );
    }
    if( a_save != NULL ) {
        LAPACKE_free( a_save );
    }
    if( jpvt != NULL ) {
        LAPACKE_free( jpvt );
    }
    if( jpvt_i != NULL ) {
        LAPACKE_free( jpvt_i );
    }
    if( jpvt_save != NULL ) {
        LAPACKE_free( jpvt_save );
    }
    if( tau != NULL ) {
        LAPACKE_free( tau );
    }
    if( tau_i != NULL ) {
        LAPACKE_free( tau_i );
    }
    if( tau_save != NULL ) {
        LAPACKE_free( tau_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: dgeqpf scalar parameters initialization */
static void init_scalars_dgeqpf( lapack_int *m, lapack_int *n, lapack_int *lda )
{
    *m = 6;
    *n = 5;
    *lda = 8;

    return;
}

/* Auxiliary functions: dgeqpf array parameters initialization */
static void init_a( lapack_int size, double *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = 0;
    }
    a[0] = -8.99999999999999970e-002;  /* a[0,0] */
    a[8] = 1.40000000000000010e-001;  /* a[0,1] */
    a[16] = -4.60000000000000020e-001;  /* a[0,2] */
    a[24] = 6.80000000000000050e-001;  /* a[0,3] */
    a[32] = 1.29000000000000000e+000;  /* a[0,4] */
    a[1] = -1.56000000000000010e+000;  /* a[1,0] */
    a[9] = 2.00000000000000010e-001;  /* a[1,1] */
    a[17] = 2.89999999999999980e-001;  /* a[1,2] */
    a[25] = 1.09000000000000010e+000;  /* a[1,3] */
    a[33] = 5.10000000000000010e-001;  /* a[1,4] */
    a[2] = -1.48000000000000000e+000;  /* a[2,0] */
    a[10] = -4.29999999999999990e-001;  /* a[2,1] */
    a[18] = 8.90000000000000010e-001;  /* a[2,2] */
    a[26] = -7.09999999999999960e-001;  /* a[2,3] */
    a[34] = -9.59999999999999960e-001;  /* a[2,4] */
    a[3] = -1.09000000000000010e+000;  /* a[3,0] */
    a[11] = 8.39999999999999970e-001;  /* a[3,1] */
    a[19] = 7.70000000000000020e-001;  /* a[3,2] */
    a[27] = 2.10999999999999990e+000;  /* a[3,3] */
    a[35] = -1.27000000000000000e+000;  /* a[3,4] */
    a[4] = 8.00000000000000020e-002;  /* a[4,0] */
    a[12] = 5.50000000000000040e-001;  /* a[4,1] */
    a[20] = -1.12999999999999990e+000;  /* a[4,2] */
    a[28] = 1.40000000000000010e-001;  /* a[4,3] */
    a[36] = 1.74000000000000000e+000;  /* a[4,4] */
    a[5] = -1.59000000000000010e+000;  /* a[5,0] */
    a[13] = -7.19999999999999970e-001;  /* a[5,1] */
    a[21] = 1.06000000000000010e+000;  /* a[5,2] */
    a[29] = 1.24000000000000000e+000;  /* a[5,3] */
    a[37] = 3.40000000000000020e-001;  /* a[5,4] */
}
static void init_jpvt( lapack_int size, lapack_int *jpvt ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        jpvt[i] = 0;
    }
    jpvt[0] = 0;
    jpvt[1] = 0;
    jpvt[2] = 0;
    jpvt[3] = 0;
    jpvt[4] = 0;
}
static void init_tau( lapack_int size, double *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = 0;
    }
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to dgeqpf results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dgeqpf( double *a, double *a_i, lapack_int *jpvt,
                           lapack_int *jpvt_i, double *tau, double *tau_i,
                           lapack_int info, lapack_int info_i, lapack_int lda,
                           lapack_int m, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < lda*n; i++ ) {
        failed += compare_doubles(a[i],a_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += (jpvt[i] == jpvt_i[i]) ? 0 : 1;
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        failed += compare_doubles(tau[i],tau_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
