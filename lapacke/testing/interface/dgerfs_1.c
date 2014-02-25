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
* dgerfs_1 is the test program for the C interface to LAPACK
* routine dgerfs
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

static void init_scalars_dgerfs( char *trans, lapack_int *n, lapack_int *nrhs,
                                 lapack_int *lda, lapack_int *ldaf,
                                 lapack_int *ldb, lapack_int *ldx );
static void init_a( lapack_int size, double *a );
static void init_af( lapack_int size, double *af );
static void init_ipiv( lapack_int size, lapack_int *ipiv );
static void init_b( lapack_int size, double *b );
static void init_x( lapack_int size, double *x );
static void init_ferr( lapack_int size, double *ferr );
static void init_berr( lapack_int size, double *berr );
static void init_work( lapack_int size, double *work );
static void init_iwork( lapack_int size, lapack_int *iwork );
static int compare_dgerfs( double *x, double *x_i, double *ferr, double *ferr_i,
                           double *berr, double *berr_i, lapack_int info,
                           lapack_int info_i, lapack_int ldx, lapack_int nrhs );

int main(void)
{
    /* Local scalars */
    char trans, trans_i;
    lapack_int n, n_i;
    lapack_int nrhs, nrhs_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    lapack_int ldaf, ldaf_i;
    lapack_int ldaf_r;
    lapack_int ldb, ldb_i;
    lapack_int ldb_r;
    lapack_int ldx, ldx_i;
    lapack_int ldx_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *a = NULL, *a_i = NULL;
    double *af = NULL, *af_i = NULL;
    lapack_int *ipiv = NULL, *ipiv_i = NULL;
    double *b = NULL, *b_i = NULL;
    double *x = NULL, *x_i = NULL;
    double *ferr = NULL, *ferr_i = NULL;
    double *berr = NULL, *berr_i = NULL;
    double *work = NULL, *work_i = NULL;
    lapack_int *iwork = NULL, *iwork_i = NULL;
    double *x_save = NULL;
    double *ferr_save = NULL;
    double *berr_save = NULL;
    double *a_r = NULL;
    double *af_r = NULL;
    double *b_r = NULL;
    double *x_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dgerfs( &trans, &n, &nrhs, &lda, &ldaf, &ldb, &ldx );
    lda_r = n+2;
    ldaf_r = n+2;
    ldb_r = nrhs+2;
    ldx_r = nrhs+2;
    trans_i = trans;
    n_i = n;
    nrhs_i = nrhs;
    lda_i = lda;
    ldaf_i = ldaf;
    ldb_i = ldb;
    ldx_i = ldx;

    /* Allocate memory for the LAPACK routine arrays */
    a = (double *)LAPACKE_malloc( lda*n * sizeof(double) );
    af = (double *)LAPACKE_malloc( ldaf*n * sizeof(double) );
    ipiv = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b = (double *)LAPACKE_malloc( ldb*nrhs * sizeof(double) );
    x = (double *)LAPACKE_malloc( ldx*nrhs * sizeof(double) );
    ferr = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    berr = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    work = (double *)LAPACKE_malloc( 3*n * sizeof(double) );
    iwork = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    a_i = (double *)LAPACKE_malloc( lda*n * sizeof(double) );
    af_i = (double *)LAPACKE_malloc( ldaf*n * sizeof(double) );
    ipiv_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b_i = (double *)LAPACKE_malloc( ldb*nrhs * sizeof(double) );
    x_i = (double *)LAPACKE_malloc( ldx*nrhs * sizeof(double) );
    ferr_i = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    berr_i = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( 3*n * sizeof(double) );
    iwork_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    x_save = (double *)LAPACKE_malloc( ldx*nrhs * sizeof(double) );
    ferr_save = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    berr_save = (double *)LAPACKE_malloc( nrhs * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    a_r = (double *)LAPACKE_malloc( n*(n+2) * sizeof(double) );
    af_r = (double *)LAPACKE_malloc( n*(n+2) * sizeof(double) );
    b_r = (double *)LAPACKE_malloc( n*(nrhs+2) * sizeof(double) );
    x_r = (double *)LAPACKE_malloc( n*(nrhs+2) * sizeof(double) );

    /* Initialize input arrays */
    init_a( lda*n, a );
    init_af( ldaf*n, af );
    init_ipiv( n, ipiv );
    init_b( ldb*nrhs, b );
    init_x( ldx*nrhs, x );
    init_ferr( nrhs, ferr );
    init_berr( nrhs, berr );
    init_work( 3*n, work );
    init_iwork( n, iwork );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_save[i] = x[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        ferr_save[i] = ferr[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        berr_save[i] = berr[i];
    }

    /* Call the LAPACK routine */
    dgerfs_( &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx,
             ferr, berr, work, iwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < ldaf*n; i++ ) {
        af_i[i] = af[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_i[i] = x_save[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        ferr_i[i] = ferr_save[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        berr_i[i] = berr_save[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_dgerfs_work( LAPACK_COL_MAJOR, trans_i, n_i, nrhs_i, a_i,
                                  lda_i, af_i, ldaf_i, ipiv_i, b_i, ldb_i, x_i,
                                  ldx_i, ferr_i, berr_i, work_i, iwork_i );

    failed = compare_dgerfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dgerfs\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dgerfs\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < ldaf*n; i++ ) {
        af_i[i] = af[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_i[i] = x_save[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        ferr_i[i] = ferr_save[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        berr_i[i] = berr_save[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_dgerfs( LAPACK_COL_MAJOR, trans_i, n_i, nrhs_i, a_i, lda_i,
                             af_i, ldaf_i, ipiv_i, b_i, ldb_i, x_i, ldx_i,
                             ferr_i, berr_i );

    failed = compare_dgerfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dgerfs\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dgerfs\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < ldaf*n; i++ ) {
        af_i[i] = af[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_i[i] = x_save[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        ferr_i[i] = ferr_save[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        berr_i[i] = berr_save[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, af_i, ldaf, af_r, n+2 );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_dgerfs_work( LAPACK_ROW_MAJOR, trans_i, n_i, nrhs_i, a_r,
                                  lda_r, af_r, ldaf_r, ipiv_i, b_r, ldb_r, x_r,
                                  ldx_r, ferr_i, berr_i, work_i, iwork_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_dgerfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dgerfs\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dgerfs\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < ldaf*n; i++ ) {
        af_i[i] = af[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_i[i] = x_save[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        ferr_i[i] = ferr_save[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        berr_i[i] = berr_save[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, af_i, ldaf, af_r, n+2 );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_dgerfs( LAPACK_ROW_MAJOR, trans_i, n_i, nrhs_i, a_r, lda_r,
                             af_r, ldaf_r, ipiv_i, b_r, ldb_r, x_r, ldx_r,
                             ferr_i, berr_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_dgerfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dgerfs\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dgerfs\n" );
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
    if( af != NULL ) {
        LAPACKE_free( af );
    }
    if( af_i != NULL ) {
        LAPACKE_free( af_i );
    }
    if( af_r != NULL ) {
        LAPACKE_free( af_r );
    }
    if( ipiv != NULL ) {
        LAPACKE_free( ipiv );
    }
    if( ipiv_i != NULL ) {
        LAPACKE_free( ipiv_i );
    }
    if( b != NULL ) {
        LAPACKE_free( b );
    }
    if( b_i != NULL ) {
        LAPACKE_free( b_i );
    }
    if( b_r != NULL ) {
        LAPACKE_free( b_r );
    }
    if( x != NULL ) {
        LAPACKE_free( x );
    }
    if( x_i != NULL ) {
        LAPACKE_free( x_i );
    }
    if( x_r != NULL ) {
        LAPACKE_free( x_r );
    }
    if( x_save != NULL ) {
        LAPACKE_free( x_save );
    }
    if( ferr != NULL ) {
        LAPACKE_free( ferr );
    }
    if( ferr_i != NULL ) {
        LAPACKE_free( ferr_i );
    }
    if( ferr_save != NULL ) {
        LAPACKE_free( ferr_save );
    }
    if( berr != NULL ) {
        LAPACKE_free( berr );
    }
    if( berr_i != NULL ) {
        LAPACKE_free( berr_i );
    }
    if( berr_save != NULL ) {
        LAPACKE_free( berr_save );
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

/* Auxiliary function: dgerfs scalar parameters initialization */
static void init_scalars_dgerfs( char *trans, lapack_int *n, lapack_int *nrhs,
                                 lapack_int *lda, lapack_int *ldaf,
                                 lapack_int *ldb, lapack_int *ldx )
{
    *trans = 'N';
    *n = 4;
    *nrhs = 2;
    *lda = 8;
    *ldaf = 8;
    *ldb = 8;
    *ldx = 8;

    return;
}

/* Auxiliary functions: dgerfs array parameters initialization */
static void init_a( lapack_int size, double *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = 0;
    }
    a[0] = 1.80000000000000000e+000;  /* a[0,0] */
    a[8] = 2.87999999999999990e+000;  /* a[0,1] */
    a[16] = 2.04999999999999980e+000;  /* a[0,2] */
    a[24] = -8.90000000000000010e-001;  /* a[0,3] */
    a[1] = 5.25000000000000000e+000;  /* a[1,0] */
    a[9] = -2.95000000000000020e+000;  /* a[1,1] */
    a[17] = -9.49999999999999960e-001;  /* a[1,2] */
    a[25] = -3.79999999999999980e+000;  /* a[1,3] */
    a[2] = 1.58000000000000010e+000;  /* a[2,0] */
    a[10] = -2.68999999999999990e+000;  /* a[2,1] */
    a[18] = -2.89999999999999990e+000;  /* a[2,2] */
    a[26] = -1.04000000000000000e+000;  /* a[2,3] */
    a[3] = -1.11000000000000010e+000;  /* a[3,0] */
    a[11] = -6.60000000000000030e-001;  /* a[3,1] */
    a[19] = -5.89999999999999970e-001;  /* a[3,2] */
    a[27] = 8.00000000000000040e-001;  /* a[3,3] */
}
static void init_af( lapack_int size, double *af ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        af[i] = 0;
    }
    af[0] = 5.25000000000000000e+000;  /* af[0,0] */
    af[8] = -2.95000000000000020e+000;  /* af[0,1] */
    af[16] = -9.49999999999999960e-001;  /* af[0,2] */
    af[24] = -3.79999999999999980e+000;  /* af[0,3] */
    af[1] = 3.42857142857142860e-001;  /* af[1,0] */
    af[9] = 3.89142857142857150e+000;  /* af[1,1] */
    af[17] = 2.37571428571428540e+000;  /* af[1,2] */
    af[25] = 4.12857142857142700e-001;  /* af[1,3] */
    af[2] = 3.00952380952380970e-001;  /* af[2,0] */
    af[10] = -4.63117963778756640e-001;  /* af[2,1] */
    af[18] = -1.51385927557513480e+000;  /* af[2,2] */
    af[26] = 2.94820606950562780e-001;  /* af[2,3] */
    af[3] = -2.11428571428571440e-001;  /* af[3,0] */
    af[11] = -3.29882525697503700e-001;  /* af[3,1] */
    af[19] = 4.72336766398369880e-003;  /* af[3,2] */
    af[27] = 1.31373239487851680e-001;  /* af[3,3] */
}
static void init_ipiv( lapack_int size, lapack_int *ipiv ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ipiv[i] = 0;
    }
    ipiv[0] = 2;
    ipiv[1] = 2;
    ipiv[2] = 3;
    ipiv[3] = 4;
}
static void init_b( lapack_int size, double *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = 0;
    }
    b[0] = 9.51999999999999960e+000;  /* b[0,0] */
    b[8] = 1.84699999999999990e+001;  /* b[0,1] */
    b[1] = 2.43500000000000010e+001;  /* b[1,0] */
    b[9] = 2.25000000000000000e+000;  /* b[1,1] */
    b[2] = 7.70000000000000020e-001;  /* b[2,0] */
    b[10] = -1.32799999999999990e+001;  /* b[2,1] */
    b[3] = -6.21999999999999980e+000;  /* b[3,0] */
    b[11] = -6.21000000000000000e+000;  /* b[3,1] */
}
static void init_x( lapack_int size, double *x ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        x[i] = 0;
    }
    x[0] = 1.00000000000000270e+000;  /* x[0,0] */
    x[8] = 3.00000000000000000e+000;  /* x[0,1] */
    x[1] = -1.00000000000000200e+000;  /* x[1,0] */
    x[9] = 1.99999999999999960e+000;  /* x[1,1] */
    x[2] = 3.00000000000000220e+000;  /* x[2,0] */
    x[10] = 4.00000000000000090e+000;  /* x[2,1] */
    x[3] = -4.99999999999999560e+000;  /* x[3,0] */
    x[11] = 1.00000000000000000e+000;  /* x[3,1] */
}
static void init_ferr( lapack_int size, double *ferr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ferr[i] = 0;
    }
}
static void init_berr( lapack_int size, double *berr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        berr[i] = 0;
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

/* Auxiliary function: C interface to dgerfs results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dgerfs( double *x, double *x_i, double *ferr, double *ferr_i,
                           double *berr, double *berr_i, lapack_int info,
                           lapack_int info_i, lapack_int ldx, lapack_int nrhs )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldx*nrhs; i++ ) {
        failed += compare_doubles(x[i],x_i[i]);
    }
    for( i = 0; i < nrhs; i++ ) {
        failed += compare_doubles(ferr[i],ferr_i[i]);
    }
    for( i = 0; i < nrhs; i++ ) {
        failed += compare_doubles(berr[i],berr_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
