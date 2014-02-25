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
* cgbrfs_1 is the test program for the C interface to LAPACK
* routine cgbrfs
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

static void init_scalars_cgbrfs( char *trans, lapack_int *n, lapack_int *kl,
                                 lapack_int *ku, lapack_int *nrhs,
                                 lapack_int *ldab, lapack_int *ldafb,
                                 lapack_int *ldb, lapack_int *ldx );
static void init_ab( lapack_int size, lapack_complex_float *ab );
static void init_afb( lapack_int size, lapack_complex_float *afb );
static void init_ipiv( lapack_int size, lapack_int *ipiv );
static void init_b( lapack_int size, lapack_complex_float *b );
static void init_x( lapack_int size, lapack_complex_float *x );
static void init_ferr( lapack_int size, float *ferr );
static void init_berr( lapack_int size, float *berr );
static void init_work( lapack_int size, lapack_complex_float *work );
static void init_rwork( lapack_int size, float *rwork );
static int compare_cgbrfs( lapack_complex_float *x, lapack_complex_float *x_i,
                           float *ferr, float *ferr_i, float *berr,
                           float *berr_i, lapack_int info, lapack_int info_i,
                           lapack_int ldx, lapack_int nrhs );

int main(void)
{
    /* Local scalars */
    char trans, trans_i;
    lapack_int n, n_i;
    lapack_int kl, kl_i;
    lapack_int ku, ku_i;
    lapack_int nrhs, nrhs_i;
    lapack_int ldab, ldab_i;
    lapack_int ldab_r;
    lapack_int ldafb, ldafb_i;
    lapack_int ldafb_r;
    lapack_int ldb, ldb_i;
    lapack_int ldb_r;
    lapack_int ldx, ldx_i;
    lapack_int ldx_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *ab = NULL, *ab_i = NULL;
    lapack_complex_float *afb = NULL, *afb_i = NULL;
    lapack_int *ipiv = NULL, *ipiv_i = NULL;
    lapack_complex_float *b = NULL, *b_i = NULL;
    lapack_complex_float *x = NULL, *x_i = NULL;
    float *ferr = NULL, *ferr_i = NULL;
    float *berr = NULL, *berr_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    float *rwork = NULL, *rwork_i = NULL;
    lapack_complex_float *x_save = NULL;
    float *ferr_save = NULL;
    float *berr_save = NULL;
    lapack_complex_float *ab_r = NULL;
    lapack_complex_float *afb_r = NULL;
    lapack_complex_float *b_r = NULL;
    lapack_complex_float *x_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_cgbrfs( &trans, &n, &kl, &ku, &nrhs, &ldab, &ldafb, &ldb,
                         &ldx );
    ldab_r = n+2;
    ldafb_r = n+2;
    ldb_r = nrhs+2;
    ldx_r = nrhs+2;
    trans_i = trans;
    n_i = n;
    kl_i = kl;
    ku_i = ku;
    nrhs_i = nrhs;
    ldab_i = ldab;
    ldafb_i = ldafb;
    ldb_i = ldb;
    ldx_i = ldx;

    /* Allocate memory for the LAPACK routine arrays */
    ab = (lapack_complex_float *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_float) );
    afb = (lapack_complex_float *)
        LAPACKE_malloc( ldafb*n * sizeof(lapack_complex_float) );
    ipiv = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b = (lapack_complex_float *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_float) );
    x = (lapack_complex_float *)
        LAPACKE_malloc( ldx*nrhs * sizeof(lapack_complex_float) );
    ferr = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    berr = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_float) );
    rwork = (float *)LAPACKE_malloc( n * sizeof(float) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (lapack_complex_float *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_float) );
    afb_i = (lapack_complex_float *)
        LAPACKE_malloc( ldafb*n * sizeof(lapack_complex_float) );
    ipiv_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b_i = (lapack_complex_float *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_float) );
    x_i = (lapack_complex_float *)
        LAPACKE_malloc( ldx*nrhs * sizeof(lapack_complex_float) );
    ferr_i = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    berr_i = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_float) );
    rwork_i = (float *)LAPACKE_malloc( n * sizeof(float) );

    /* Allocate memory for the backup arrays */
    x_save = (lapack_complex_float *)
        LAPACKE_malloc( ldx*nrhs * sizeof(lapack_complex_float) );
    ferr_save = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    berr_save = (float *)LAPACKE_malloc( nrhs * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    ab_r = (lapack_complex_float *)
        LAPACKE_malloc( (kl+ku+1)*(n+2) * sizeof(lapack_complex_float) );
    afb_r = (lapack_complex_float *)
        LAPACKE_malloc( ((2*kl+ku+1)*(n+2)) * sizeof(lapack_complex_float) );
    b_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(nrhs+2) * sizeof(lapack_complex_float) );
    x_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(nrhs+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_afb( ldafb*n, afb );
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
    cgbrfs_( &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb,
             x, &ldx, ferr, berr, work, rwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldafb*n; i++ ) {
        afb_i[i] = afb[i];
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
    info_i = LAPACKE_cgbrfs_work( LAPACK_COL_MAJOR, trans_i, n_i, kl_i, ku_i,
                                  nrhs_i, ab_i, ldab_i, afb_i, ldafb_i, ipiv_i,
                                  b_i, ldb_i, x_i, ldx_i, ferr_i, berr_i,
                                  work_i, rwork_i );

    failed = compare_cgbrfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to cgbrfs\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to cgbrfs\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldafb*n; i++ ) {
        afb_i[i] = afb[i];
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
    info_i = LAPACKE_cgbrfs( LAPACK_COL_MAJOR, trans_i, n_i, kl_i, ku_i, nrhs_i,
                             ab_i, ldab_i, afb_i, ldafb_i, ipiv_i, b_i, ldb_i,
                             x_i, ldx_i, ferr_i, berr_i );

    failed = compare_cgbrfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to cgbrfs\n" );
    } else {
        printf( "FAILED: column-major high-level interface to cgbrfs\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldafb*n; i++ ) {
        afb_i[i] = afb[i];
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

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, 2*kl+ku+1, n, afb_i, ldafb, afb_r,
                       n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_cgbrfs_work( LAPACK_ROW_MAJOR, trans_i, n_i, kl_i, ku_i,
                                  nrhs_i, ab_r, ldab_r, afb_r, ldafb_r, ipiv_i,
                                  b_r, ldb_r, x_r, ldx_r, ferr_i, berr_i,
                                  work_i, rwork_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_cgbrfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to cgbrfs\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to cgbrfs\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldafb*n; i++ ) {
        afb_i[i] = afb[i];
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
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, 2*kl+ku+1, n, afb_i, ldafb, afb_r,
                       n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_cgbrfs( LAPACK_ROW_MAJOR, trans_i, n_i, kl_i, ku_i, nrhs_i,
                             ab_r, ldab_r, afb_r, ldafb_r, ipiv_i, b_r, ldb_r,
                             x_r, ldx_r, ferr_i, berr_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_cgbrfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to cgbrfs\n" );
    } else {
        printf( "FAILED: row-major high-level interface to cgbrfs\n" );
    }

    /* Release memory */
    if( ab != NULL ) {
        LAPACKE_free( ab );
    }
    if( ab_i != NULL ) {
        LAPACKE_free( ab_i );
    }
    if( ab_r != NULL ) {
        LAPACKE_free( ab_r );
    }
    if( afb != NULL ) {
        LAPACKE_free( afb );
    }
    if( afb_i != NULL ) {
        LAPACKE_free( afb_i );
    }
    if( afb_r != NULL ) {
        LAPACKE_free( afb_r );
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

/* Auxiliary function: cgbrfs scalar parameters initialization */
static void init_scalars_cgbrfs( char *trans, lapack_int *n, lapack_int *kl,
                                 lapack_int *ku, lapack_int *nrhs,
                                 lapack_int *ldab, lapack_int *ldafb,
                                 lapack_int *ldb, lapack_int *ldx )
{
    *trans = 'N';
    *n = 4;
    *kl = 1;
    *ku = 2;
    *nrhs = 2;
    *ldab = 17;
    *ldafb = 25;
    *ldb = 8;
    *ldx = 8;

    return;
}

/* Auxiliary functions: cgbrfs array parameters initialization */
static void init_ab( lapack_int size, lapack_complex_float *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    ab[0] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[17] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[34] = lapack_make_complex_float( 9.700000286e-001, -2.839999914e+000 );
    ab[51] = lapack_make_complex_float( 5.899999738e-001, -4.799999893e-001 );
    ab[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[18] = lapack_make_complex_float( -2.049999952e+000, -8.500000238e-001 );
    ab[35] = lapack_make_complex_float( -3.990000010e+000, 4.010000229e+000 );
    ab[52] = lapack_make_complex_float( 3.329999924e+000, -1.039999962e+000 );
    ab[2] = lapack_make_complex_float( -1.649999976e+000, 2.259999990e+000 );
    ab[19] = lapack_make_complex_float( -1.480000019e+000, -1.750000000e+000 );
    ab[36] = lapack_make_complex_float( -1.059999943e+000, 1.940000057e+000 );
    ab[53] = lapack_make_complex_float( -4.600000083e-001, -1.720000029e+000 );
    ab[3] = lapack_make_complex_float( 0.000000000e+000, 6.300000191e+000 );
    ab[20] = lapack_make_complex_float( -7.699999809e-001, 2.829999924e+000 );
    ab[37] = lapack_make_complex_float( 4.480000019e+000, -1.090000033e+000 );
    ab[54] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_afb( lapack_int size, lapack_complex_float *afb ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        afb[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    afb[0] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    afb[25] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    afb[50] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    afb[75] = lapack_make_complex_float( 5.899999738e-001, -4.799999893e-001 );
    afb[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    afb[26] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    afb[51] = lapack_make_complex_float( -3.990000010e+000, 4.010000229e+000 );
    afb[76] = lapack_make_complex_float( 3.329999924e+000, -1.039999962e+000 );
    afb[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    afb[27] = lapack_make_complex_float( -1.480000019e+000, -1.750000000e+000 );
    afb[52] = lapack_make_complex_float( -1.059999943e+000, 1.940000057e+000 );
    afb[77] = lapack_make_complex_float( -1.769209266e+000, -1.858747363e+000 );
    afb[3] = lapack_make_complex_float( 0.000000000e+000, 6.300000191e+000 );
    afb[28] = lapack_make_complex_float( -7.699999809e-001, 2.829999924e+000 );
    afb[53] = lapack_make_complex_float( 4.930266857e+000, -3.008563757e+000 );
    afb[78] = lapack_make_complex_float( 4.337749183e-001, 1.232528687e-001 );
    afb[4] = lapack_make_complex_float( 3.587301373e-001, 2.619047463e-001 );
    afb[29] = lapack_make_complex_float( 2.314260751e-001, 6.357648969e-001 );
    afb[54] = lapack_make_complex_float( 7.604227066e-001, 2.429442555e-001 );
    afb[79] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_ipiv( lapack_int size, lapack_int *ipiv ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ipiv[i] = 0;
    }
    ipiv[0] = 2;
    ipiv[1] = 3;
    ipiv[2] = 3;
    ipiv[3] = 4;
}
static void init_b( lapack_int size, lapack_complex_float *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    b[0] = lapack_make_complex_float( -1.059999943e+000, 2.150000000e+001 );
    b[8] = lapack_make_complex_float( 1.285000038e+001, 2.839999914e+000 );
    b[1] = lapack_make_complex_float( -2.271999931e+001, -5.390000153e+001 );
    b[9] = lapack_make_complex_float( -7.022000122e+001, 2.156999969e+001 );
    b[2] = lapack_make_complex_float( 2.823999977e+001, -3.859999847e+001 );
    b[10] = lapack_make_complex_float( -2.072999954e+001, -1.230000019e+000 );
    b[3] = lapack_make_complex_float( -3.456000137e+001, 1.672999954e+001 );
    b[11] = lapack_make_complex_float( 2.601000023e+001, 3.196999931e+001 );
}
static void init_x( lapack_int size, lapack_complex_float *x ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        x[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    x[0] = lapack_make_complex_float( -3.000000715e+000, 1.999998450e+000 );
    x[8] = lapack_make_complex_float( 9.999980330e-001, 5.999998093e+000 );
    x[1] = lapack_make_complex_float( 1.000001788e+000, -7.000001431e+000 );
    x[9] = lapack_make_complex_float( -6.999998093e+000, -4.000003815e+000 );
    x[2] = lapack_make_complex_float( -4.999999523e+000, 4.000000000e+000 );
    x[10] = lapack_make_complex_float( 3.000001431e+000, 4.999998569e+000 );
    x[3] = lapack_make_complex_float( 6.000000000e+000, -8.000000954e+000 );
    x[11] = lapack_make_complex_float( -8.000001907e+000, 1.999996185e+000 );
}
static void init_ferr( lapack_int size, float *ferr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ferr[i] = 0;
    }
}
static void init_berr( lapack_int size, float *berr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        berr[i] = 0;
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

/* Auxiliary function: C interface to cgbrfs results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_cgbrfs( lapack_complex_float *x, lapack_complex_float *x_i,
                           float *ferr, float *ferr_i, float *berr,
                           float *berr_i, lapack_int info, lapack_int info_i,
                           lapack_int ldx, lapack_int nrhs )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldx*nrhs; i++ ) {
        failed += compare_complex_floats(x[i],x_i[i]);
    }
    for( i = 0; i < nrhs; i++ ) {
        failed += compare_floats(ferr[i],ferr_i[i]);
    }
    for( i = 0; i < nrhs; i++ ) {
        failed += compare_floats(berr[i],berr_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
