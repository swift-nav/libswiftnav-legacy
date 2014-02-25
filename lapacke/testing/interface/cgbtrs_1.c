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
* cgbtrs_1 is the test program for the C interface to LAPACK
* routine cgbtrs
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

static void init_scalars_cgbtrs( char *trans, lapack_int *n, lapack_int *kl,
                                 lapack_int *ku, lapack_int *nrhs,
                                 lapack_int *ldab, lapack_int *ldb );
static void init_ab( lapack_int size, lapack_complex_float *ab );
static void init_ipiv( lapack_int size, lapack_int *ipiv );
static void init_b( lapack_int size, lapack_complex_float *b );
static int compare_cgbtrs( lapack_complex_float *b, lapack_complex_float *b_i,
                           lapack_int info, lapack_int info_i, lapack_int ldb,
                           lapack_int nrhs );

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
    lapack_int ldb, ldb_i;
    lapack_int ldb_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *ab = NULL, *ab_i = NULL;
    lapack_int *ipiv = NULL, *ipiv_i = NULL;
    lapack_complex_float *b = NULL, *b_i = NULL;
    lapack_complex_float *b_save = NULL;
    lapack_complex_float *ab_r = NULL;
    lapack_complex_float *b_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_cgbtrs( &trans, &n, &kl, &ku, &nrhs, &ldab, &ldb );
    ldab_r = n+2;
    ldb_r = nrhs+2;
    trans_i = trans;
    n_i = n;
    kl_i = kl;
    ku_i = ku;
    nrhs_i = nrhs;
    ldab_i = ldab;
    ldb_i = ldb;

    /* Allocate memory for the LAPACK routine arrays */
    ab = (lapack_complex_float *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_float) );
    ipiv = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b = (lapack_complex_float *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (lapack_complex_float *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_float) );
    ipiv_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b_i = (lapack_complex_float *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    b_save = (lapack_complex_float *)
        LAPACKE_malloc( ldb*nrhs * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    ab_r = (lapack_complex_float *)
        LAPACKE_malloc( ((2*kl+ku+1)*(n+2)) * sizeof(lapack_complex_float) );
    b_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(nrhs+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_ipiv( n, ipiv );
    init_b( ldb*nrhs, b );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_save[i] = b[i];
    }

    /* Call the LAPACK routine */
    cgbtrs_( &trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b_save[i];
    }
    info_i = LAPACKE_cgbtrs_work( LAPACK_COL_MAJOR, trans_i, n_i, kl_i, ku_i,
                                  nrhs_i, ab_i, ldab_i, ipiv_i, b_i, ldb_i );

    failed = compare_cgbtrs( b, b_i, info, info_i, ldb, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to cgbtrs\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to cgbtrs\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b_save[i];
    }
    info_i = LAPACKE_cgbtrs( LAPACK_COL_MAJOR, trans_i, n_i, kl_i, ku_i, nrhs_i,
                             ab_i, ldab_i, ipiv_i, b_i, ldb_i );

    failed = compare_cgbtrs( b, b_i, info, info_i, ldb, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to cgbtrs\n" );
    } else {
        printf( "FAILED: column-major high-level interface to cgbtrs\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b_save[i];
    }

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, 2*kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    info_i = LAPACKE_cgbtrs_work( LAPACK_ROW_MAJOR, trans_i, n_i, kl_i, ku_i,
                                  nrhs_i, ab_r, ldab_r, ipiv_i, b_r, ldb_r );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, nrhs, b_r, nrhs+2, b_i, ldb );

    failed = compare_cgbtrs( b, b_i, info, info_i, ldb, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to cgbtrs\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to cgbtrs\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, 2*kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    info_i = LAPACKE_cgbtrs( LAPACK_ROW_MAJOR, trans_i, n_i, kl_i, ku_i, nrhs_i,
                             ab_r, ldab_r, ipiv_i, b_r, ldb_r );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, nrhs, b_r, nrhs+2, b_i, ldb );

    failed = compare_cgbtrs( b, b_i, info, info_i, ldb, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to cgbtrs\n" );
    } else {
        printf( "FAILED: row-major high-level interface to cgbtrs\n" );
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
    if( b_save != NULL ) {
        LAPACKE_free( b_save );
    }

    return 0;
}

/* Auxiliary function: cgbtrs scalar parameters initialization */
static void init_scalars_cgbtrs( char *trans, lapack_int *n, lapack_int *kl,
                                 lapack_int *ku, lapack_int *nrhs,
                                 lapack_int *ldab, lapack_int *ldb )
{
    *trans = 'N';
    *n = 4;
    *kl = 1;
    *ku = 2;
    *nrhs = 2;
    *ldab = 25;
    *ldb = 8;

    return;
}

/* Auxiliary functions: cgbtrs array parameters initialization */
static void init_ab( lapack_int size, lapack_complex_float *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    ab[0] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[25] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[50] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[75] = lapack_make_complex_float( 5.899999738e-001, -4.799999893e-001 );
    ab[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[26] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[51] = lapack_make_complex_float( -3.990000010e+000, 4.010000229e+000 );
    ab[76] = lapack_make_complex_float( 3.329999924e+000, -1.039999962e+000 );
    ab[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[27] = lapack_make_complex_float( -1.480000019e+000, -1.750000000e+000 );
    ab[52] = lapack_make_complex_float( -1.059999943e+000, 1.940000057e+000 );
    ab[77] = lapack_make_complex_float( -1.769209266e+000, -1.858747363e+000 );
    ab[3] = lapack_make_complex_float( 0.000000000e+000, 6.300000191e+000 );
    ab[28] = lapack_make_complex_float( -7.699999809e-001, 2.829999924e+000 );
    ab[53] = lapack_make_complex_float( 4.930266857e+000, -3.008563757e+000 );
    ab[78] = lapack_make_complex_float( 4.337749183e-001, 1.232528687e-001 );
    ab[4] = lapack_make_complex_float( 3.587301373e-001, 2.619047463e-001 );
    ab[29] = lapack_make_complex_float( 2.314260751e-001, 6.357648969e-001 );
    ab[54] = lapack_make_complex_float( 7.604227066e-001, 2.429442555e-001 );
    ab[79] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
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

/* Auxiliary function: C interface to cgbtrs results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_cgbtrs( lapack_complex_float *b, lapack_complex_float *b_i,
                           lapack_int info, lapack_int info_i, lapack_int ldb,
                           lapack_int nrhs )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldb*nrhs; i++ ) {
        failed += compare_complex_floats(b[i],b_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
