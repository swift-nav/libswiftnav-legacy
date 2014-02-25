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
* dopmtr_1 is the test program for the C interface to LAPACK
* routine dopmtr
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

static void init_scalars_dopmtr( char *side, char *uplo, char *trans,
                                 lapack_int *m, lapack_int *n,
                                 lapack_int *ldc );
static void init_ap( lapack_int size, double *ap );
static void init_tau( lapack_int size, double *tau );
static void init_c( lapack_int size, double *c );
static void init_work( lapack_int size, double *work );
static int compare_dopmtr( double *c, double *c_i, lapack_int info,
                           lapack_int info_i, lapack_int ldc, lapack_int n );

int main(void)
{
    /* Local scalars */
    char side, side_i;
    char uplo, uplo_i;
    char trans, trans_i;
    lapack_int m, m_i;
    lapack_int n, n_i;
    lapack_int ldc, ldc_i;
    lapack_int ldc_r;
    lapack_int info, info_i;
    /* Declare scalars */
    lapack_int lwork;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *ap = NULL, *ap_i = NULL;
    double *tau = NULL, *tau_i = NULL;
    double *c = NULL, *c_i = NULL;
    double *work = NULL, *work_i = NULL;
    double *c_save = NULL;
    double *ap_r = NULL;
    double *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dopmtr( &side, &uplo, &trans, &m, &n, &ldc );
    lwork = MAX(m,n);
    ldc_r = n+2;
    side_i = side;
    uplo_i = uplo;
    trans_i = trans;
    m_i = m;
    n_i = n;
    ldc_i = ldc;

    /* Allocate memory for the LAPACK routine arrays */
    ap = (double *)LAPACKE_malloc( ((m*(m+1)/2)) * sizeof(double) );
    tau = (double *)LAPACKE_malloc( (m-1) * sizeof(double) );
    c = (double *)LAPACKE_malloc( ldc*n * sizeof(double) );
    work = (double *)LAPACKE_malloc( lwork * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    ap_i = (double *)LAPACKE_malloc( ((m*(m+1)/2)) * sizeof(double) );
    tau_i = (double *)LAPACKE_malloc( (m-1) * sizeof(double) );
    c_i = (double *)LAPACKE_malloc( ldc*n * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( lwork * sizeof(double) );

    /* Allocate memory for the backup arrays */
    c_save = (double *)LAPACKE_malloc( ldc*n * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    ap_r = (double *)LAPACKE_malloc( m*(m+1)/2 * sizeof(double) );
    c_r = (double *)LAPACKE_malloc( m*(n+2) * sizeof(double) );

    /* Initialize input arrays */
    init_ap( (m*(m+1)/2), ap );
    init_tau( (m-1), tau );
    init_c( ldc*n, c );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldc*n; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    dopmtr_( &side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (m*(m+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < (m-1); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dopmtr_work( LAPACK_COL_MAJOR, side_i, uplo_i, trans_i,
                                  m_i, n_i, ap_i, tau_i, c_i, ldc_i, work_i );

    failed = compare_dopmtr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dopmtr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dopmtr\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (m*(m+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < (m-1); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dopmtr( LAPACK_COL_MAJOR, side_i, uplo_i, trans_i, m_i,
                             n_i, ap_i, tau_i, c_i, ldc_i );

    failed = compare_dopmtr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dopmtr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dopmtr\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (m*(m+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < (m-1); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, m, ap_i, ap_r );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_dopmtr_work( LAPACK_ROW_MAJOR, side_i, uplo_i, trans_i,
                                  m_i, n_i, ap_r, tau_i, c_r, ldc_r, work_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_dopmtr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dopmtr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dopmtr\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (m*(m+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < (m-1); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, m, ap_i, ap_r );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_dopmtr( LAPACK_ROW_MAJOR, side_i, uplo_i, trans_i, m_i,
                             n_i, ap_r, tau_i, c_r, ldc_r );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_dopmtr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dopmtr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dopmtr\n" );
    }

    /* Release memory */
    if( ap != NULL ) {
        LAPACKE_free( ap );
    }
    if( ap_i != NULL ) {
        LAPACKE_free( ap_i );
    }
    if( ap_r != NULL ) {
        LAPACKE_free( ap_r );
    }
    if( tau != NULL ) {
        LAPACKE_free( tau );
    }
    if( tau_i != NULL ) {
        LAPACKE_free( tau_i );
    }
    if( c != NULL ) {
        LAPACKE_free( c );
    }
    if( c_i != NULL ) {
        LAPACKE_free( c_i );
    }
    if( c_r != NULL ) {
        LAPACKE_free( c_r );
    }
    if( c_save != NULL ) {
        LAPACKE_free( c_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: dopmtr scalar parameters initialization */
static void init_scalars_dopmtr( char *side, char *uplo, char *trans,
                                 lapack_int *m, lapack_int *n, lapack_int *ldc )
{
    *side = 'L';
    *uplo = 'L';
    *trans = 'N';
    *m = 4;
    *n = 2;
    *ldc = 8;

    return;
}

/* Auxiliary functions: dopmtr array parameters initialization */
static void init_ap( lapack_int size, double *ap ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ap[i] = 0;
    }
    ap[0] = 2.06999999999999980e+000;
    ap[1] = -5.82575317019181590e+000;
    ap[2] = 4.33179344221786720e-001;
    ap[3] = -1.18608629965489210e-001;
    ap[4] = 1.47409370819755310e+000;
    ap[5] = 2.62404517879558740e+000;
    ap[6] = 8.06288153277579080e-001;
    ap[7] = -6.49159507545784330e-001;
    ap[8] = 9.16272756321918620e-001;
    ap[9] = -1.69493420065176800e+000;
}
static void init_tau( lapack_int size, double *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = 0;
    }
    tau[0] = 1.66429178973824920e+000;
    tau[1] = 1.21204732416214210e+000;
    tau[2] = 0.00000000000000000e+000;
}
static void init_c( lapack_int size, double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = 0;
    }
    c[0] = 5.65759178822387240e-001;  /* c[0,0] */
    c[8] = -2.32842430803157250e-001;  /* c[0,1] */
    c[1] = 6.86917957250591790e-001;  /* c[1,0] */
    c[9] = -1.62617096149163560e-001;  /* c[1,1] */
    c[2] = -4.39588937213165000e-001;  /* c[2,0] */
    c[10] = -3.01727334388271980e-001;  /* c[2,1] */
    c[3] = 1.21744970593008300e-001;  /* c[3,0] */
    c[11] = 9.10110267022979150e-001;  /* c[3,1] */
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to dopmtr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dopmtr( double *c, double *c_i, lapack_int info,
                           lapack_int info_i, lapack_int ldc, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldc*n; i++ ) {
        failed += compare_doubles(c[i],c_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
