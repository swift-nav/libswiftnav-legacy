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
* zherfs_1 is the test program for the C interface to LAPACK
* routine zherfs
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

static void init_scalars_zherfs( char *uplo, lapack_int *n, lapack_int *nrhs,
                                 lapack_int *lda, lapack_int *ldaf,
                                 lapack_int *ldb, lapack_int *ldx );
static void init_a( lapack_int size, lapack_complex_double *a );
static void init_af( lapack_int size, lapack_complex_double *af );
static void init_ipiv( lapack_int size, lapack_int *ipiv );
static void init_b( lapack_int size, lapack_complex_double *b );
static void init_x( lapack_int size, lapack_complex_double *x );
static void init_ferr( lapack_int size, double *ferr );
static void init_berr( lapack_int size, double *berr );
static void init_work( lapack_int size, lapack_complex_double *work );
static void init_rwork( lapack_int size, double *rwork );
static int compare_zherfs( lapack_complex_double *x, lapack_complex_double *x_i,
                           double *ferr, double *ferr_i, double *berr,
                           double *berr_i, lapack_int info, lapack_int info_i,
                           lapack_int ldx, lapack_int nrhs );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
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
    lapack_complex_double *a = NULL, *a_i = NULL;
    lapack_complex_double *af = NULL, *af_i = NULL;
    lapack_int *ipiv = NULL, *ipiv_i = NULL;
    lapack_complex_double *b = NULL, *b_i = NULL;
    lapack_complex_double *x = NULL, *x_i = NULL;
    double *ferr = NULL, *ferr_i = NULL;
    double *berr = NULL, *berr_i = NULL;
    lapack_complex_double *work = NULL, *work_i = NULL;
    double *rwork = NULL, *rwork_i = NULL;
    lapack_complex_double *x_save = NULL;
    double *ferr_save = NULL;
    double *berr_save = NULL;
    lapack_complex_double *a_r = NULL;
    lapack_complex_double *af_r = NULL;
    lapack_complex_double *b_r = NULL;
    lapack_complex_double *x_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_zherfs( &uplo, &n, &nrhs, &lda, &ldaf, &ldb, &ldx );
    lda_r = n+2;
    ldaf_r = n+2;
    ldb_r = nrhs+2;
    ldx_r = nrhs+2;
    uplo_i = uplo;
    n_i = n;
    nrhs_i = nrhs;
    lda_i = lda;
    ldaf_i = ldaf;
    ldb_i = ldb;
    ldx_i = ldx;

    /* Allocate memory for the LAPACK routine arrays */
    a = (lapack_complex_double *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_double) );
    af = (lapack_complex_double *)
        LAPACKE_malloc( ldaf*n * sizeof(lapack_complex_double) );
    ipiv = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b = (lapack_complex_double *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_double) );
    x = (lapack_complex_double *)
        LAPACKE_malloc( ldx*nrhs * sizeof(lapack_complex_double) );
    ferr = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    berr = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    work = (lapack_complex_double *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_double) );
    rwork = (double *)LAPACKE_malloc( n * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_double *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_double) );
    af_i = (lapack_complex_double *)
        LAPACKE_malloc( ldaf*n * sizeof(lapack_complex_double) );
    ipiv_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b_i = (lapack_complex_double *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_double) );
    x_i = (lapack_complex_double *)
        LAPACKE_malloc( ldx*nrhs * sizeof(lapack_complex_double) );
    ferr_i = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    berr_i = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    work_i = (lapack_complex_double *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_double) );
    rwork_i = (double *)LAPACKE_malloc( n * sizeof(double) );

    /* Allocate memory for the backup arrays */
    x_save = (lapack_complex_double *)
        LAPACKE_malloc( ldx*nrhs * sizeof(lapack_complex_double) );
    ferr_save = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    berr_save = (double *)LAPACKE_malloc( nrhs * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );
    af_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );
    b_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(nrhs+2) * sizeof(lapack_complex_double) );
    x_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(nrhs+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_a( lda*n, a );
    init_af( ldaf*n, af );
    init_ipiv( n, ipiv );
    init_b( ldb*nrhs, b );
    init_x( ldx*nrhs, x );
    init_ferr( nrhs, ferr );
    init_berr( nrhs, berr );
    init_work( 2*n, work );
    init_rwork( n, rwork );

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
    zherfs_( &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr,
             berr, work, rwork, &info );

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
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_zherfs_work( LAPACK_COL_MAJOR, uplo_i, n_i, nrhs_i, a_i,
                                  lda_i, af_i, ldaf_i, ipiv_i, b_i, ldb_i, x_i,
                                  ldx_i, ferr_i, berr_i, work_i, rwork_i );

    failed = compare_zherfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zherfs\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zherfs\n" );
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
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_zherfs( LAPACK_COL_MAJOR, uplo_i, n_i, nrhs_i, a_i, lda_i,
                             af_i, ldaf_i, ipiv_i, b_i, ldb_i, x_i, ldx_i,
                             ferr_i, berr_i );

    failed = compare_zherfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zherfs\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zherfs\n" );
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
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, af_i, ldaf, af_r, n+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_zherfs_work( LAPACK_ROW_MAJOR, uplo_i, n_i, nrhs_i, a_r,
                                  lda_r, af_r, ldaf_r, ipiv_i, b_r, ldb_r, x_r,
                                  ldx_r, ferr_i, berr_i, work_i, rwork_i );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_zherfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zherfs\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zherfs\n" );
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
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        rwork_i[i] = rwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, af_i, ldaf, af_r, n+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_zherfs( LAPACK_ROW_MAJOR, uplo_i, n_i, nrhs_i, a_r, lda_r,
                             af_r, ldaf_r, ipiv_i, b_r, ldb_r, x_r, ldx_r,
                             ferr_i, berr_i );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_zherfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zherfs\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zherfs\n" );
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
    if( rwork != NULL ) {
        LAPACKE_free( rwork );
    }
    if( rwork_i != NULL ) {
        LAPACKE_free( rwork_i );
    }

    return 0;
}

/* Auxiliary function: zherfs scalar parameters initialization */
static void init_scalars_zherfs( char *uplo, lapack_int *n, lapack_int *nrhs,
                                 lapack_int *lda, lapack_int *ldaf,
                                 lapack_int *ldb, lapack_int *ldx )
{
    *uplo = 'L';
    *n = 4;
    *nrhs = 2;
    *lda = 8;
    *ldaf = 8;
    *ldb = 8;
    *ldx = 8;

    return;
}

/* Auxiliary functions: zherfs array parameters initialization */
static void init_a( lapack_int size, lapack_complex_double *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    a[0] = lapack_make_complex_double( -1.36000000000000010e+000,
                                       0.00000000000000000e+000 );
    a[8] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    a[16] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    a[24] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    a[1] = lapack_make_complex_double( 1.58000000000000010e+000,
                                       -9.00000000000000020e-001 );
    a[9] = lapack_make_complex_double( -8.86999999999999920e+000,
                                       0.00000000000000000e+000 );
    a[17] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    a[25] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    a[2] = lapack_make_complex_double( 2.21000000000000000e+000,
                                       2.09999999999999990e-001 );
    a[10] = lapack_make_complex_double( -1.84000000000000010e+000,
                                        2.99999999999999990e-002 );
    a[18] = lapack_make_complex_double( -4.62999999999999990e+000,
                                        0.00000000000000000e+000 );
    a[26] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    a[3] = lapack_make_complex_double( 3.91000000000000010e+000,
                                       -1.50000000000000000e+000 );
    a[11] = lapack_make_complex_double( -1.78000000000000000e+000,
                                        -1.17999999999999990e+000 );
    a[19] = lapack_make_complex_double( 1.10000000000000000e-001,
                                        -1.10000000000000000e-001 );
    a[27] = lapack_make_complex_double( -1.84000000000000010e+000,
                                        0.00000000000000000e+000 );
}
static void init_af( lapack_int size, lapack_complex_double *af ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        af[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    af[0] = lapack_make_complex_double( -1.36000000000000010e+000,
                                        0.00000000000000000e+000 );
    af[8] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    af[16] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    af[24] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    af[1] = lapack_make_complex_double( 3.91000000000000010e+000,
                                        -1.50000000000000000e+000 );
    af[9] = lapack_make_complex_double( -1.84000000000000010e+000,
                                        0.00000000000000000e+000 );
    af[17] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    af[25] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    af[2] = lapack_make_complex_double( 3.10028798127124030e-001,
                                        4.33302074396270180e-002 );
    af[10] = lapack_make_complex_double( 5.63705048650877560e-001,
                                         2.85034950151971610e-001 );
    af[18] = lapack_make_complex_double( -5.41762438729157920e+000,
                                         0.00000000000000000e+000 );
    af[26] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    af[3] = lapack_make_complex_double( -1.51812020724010200e-001,
                                        3.74295842561370500e-001 );
    af[11] = lapack_make_complex_double( 3.39658279960360960e-001,
                                         3.03145181135563540e-002 );
    af[19] = lapack_make_complex_double( 2.99724464607583510e-001,
                                         1.57826837278577770e-001 );
    af[27] = lapack_make_complex_double( -7.10280989580184220e+000,
                                         0.00000000000000000e+000 );
}
static void init_ipiv( lapack_int size, lapack_int *ipiv ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ipiv[i] = 0;
    }
    ipiv[0] = -4;
    ipiv[1] = -4;
    ipiv[2] = 3;
    ipiv[3] = 4;
}
static void init_b( lapack_int size, lapack_complex_double *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    b[0] = lapack_make_complex_double( 7.79000000000000000e+000,
                                       5.48000000000000040e+000 );
    b[8] = lapack_make_complex_double( -3.53900000000000010e+001,
                                       1.80100000000000020e+001 );
    b[1] = lapack_make_complex_double( -7.70000000000000020e-001,
                                       -1.60500000000000010e+001 );
    b[9] = lapack_make_complex_double( 4.23000000000000040e+000,
                                       -7.00199999999999960e+001 );
    b[2] = lapack_make_complex_double( -9.58000000000000010e+000,
                                       3.87999999999999990e+000 );
    b[10] = lapack_make_complex_double( -2.47899999999999990e+001,
                                        -8.40000000000000040e+000 );
    b[3] = lapack_make_complex_double( 2.98000000000000000e+000,
                                       -1.01800000000000000e+001 );
    b[11] = lapack_make_complex_double( 2.86800000000000000e+001,
                                        -3.98900000000000010e+001 );
}
static void init_x( lapack_int size, lapack_complex_double *x ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        x[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    x[0] = lapack_make_complex_double( 1.00000000000000040e+000,
                                       -1.00000000000000000e+000 );
    x[8] = lapack_make_complex_double( 3.00000000000000270e+000,
                                       -4.00000000000000000e+000 );
    x[1] = lapack_make_complex_double( -1.00000000000000020e+000,
                                       2.00000000000000040e+000 );
    x[9] = lapack_make_complex_double( -1.00000000000000040e+000,
                                       5.00000000000000090e+000 );
    x[2] = lapack_make_complex_double( 3.00000000000000000e+000,
                                       -2.00000000000000040e+000 );
    x[10] = lapack_make_complex_double( 7.00000000000000000e+000,
                                        -2.00000000000000090e+000 );
    x[3] = lapack_make_complex_double( 2.00000000000000000e+000,
                                       1.00000000000000000e+000 );
    x[11] = lapack_make_complex_double( -7.99999999999999820e+000,
                                        6.00000000000000090e+000 );
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

/* Auxiliary function: C interface to zherfs results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zherfs( lapack_complex_double *x, lapack_complex_double *x_i,
                           double *ferr, double *ferr_i, double *berr,
                           double *berr_i, lapack_int info, lapack_int info_i,
                           lapack_int ldx, lapack_int nrhs )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldx*nrhs; i++ ) {
        failed += compare_complex_doubles(x[i],x_i[i]);
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
