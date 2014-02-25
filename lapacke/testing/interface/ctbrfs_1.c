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
* ctbrfs_1 is the test program for the C interface to LAPACK
* routine ctbrfs
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

static void init_scalars_ctbrfs( char *uplo, char *trans, char *diag,
                                 lapack_int *n, lapack_int *kd,
                                 lapack_int *nrhs, lapack_int *ldab,
                                 lapack_int *ldb, lapack_int *ldx );
static void init_ab( lapack_int size, lapack_complex_float *ab );
static void init_b( lapack_int size, lapack_complex_float *b );
static void init_x( lapack_int size, lapack_complex_float *x );
static void init_ferr( lapack_int size, float *ferr );
static void init_berr( lapack_int size, float *berr );
static void init_work( lapack_int size, lapack_complex_float *work );
static void init_rwork( lapack_int size, float *rwork );
static int compare_ctbrfs( float *ferr, float *ferr_i, float *berr,
                           float *berr_i, lapack_int info, lapack_int info_i,
                           lapack_int nrhs );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
    char trans, trans_i;
    char diag, diag_i;
    lapack_int n, n_i;
    lapack_int kd, kd_i;
    lapack_int nrhs, nrhs_i;
    lapack_int ldab, ldab_i;
    lapack_int ldab_r;
    lapack_int ldb, ldb_i;
    lapack_int ldb_r;
    lapack_int ldx, ldx_i;
    lapack_int ldx_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *ab = NULL, *ab_i = NULL;
    lapack_complex_float *b = NULL, *b_i = NULL;
    lapack_complex_float *x = NULL, *x_i = NULL;
    float *ferr = NULL, *ferr_i = NULL;
    float *berr = NULL, *berr_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    float *rwork = NULL, *rwork_i = NULL;
    float *ferr_save = NULL;
    float *berr_save = NULL;
    lapack_complex_float *ab_r = NULL;
    lapack_complex_float *b_r = NULL;
    lapack_complex_float *x_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ctbrfs( &uplo, &trans, &diag, &n, &kd, &nrhs, &ldab, &ldb,
                         &ldx );
    ldab_r = n+2;
    ldb_r = nrhs+2;
    ldx_r = nrhs+2;
    uplo_i = uplo;
    trans_i = trans;
    diag_i = diag;
    n_i = n;
    kd_i = kd;
    nrhs_i = nrhs;
    ldab_i = ldab;
    ldb_i = ldb;
    ldx_i = ldx;

    /* Allocate memory for the LAPACK routine arrays */
    ab = (lapack_complex_float *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_float) );
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
    ferr_save = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    berr_save = (float *)LAPACKE_malloc( nrhs * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    ab_r = (lapack_complex_float *)
        LAPACKE_malloc( (kd+1)*(n+2) * sizeof(lapack_complex_float) );
    b_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(nrhs+2) * sizeof(lapack_complex_float) );
    x_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(nrhs+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_b( ldb*nrhs, b );
    init_x( ldx*nrhs, x );
    init_ferr( nrhs, ferr );
    init_berr( nrhs, berr );
    init_work( 2*n, work );
    init_rwork( n, rwork );

    /* Backup the ouptut arrays */
    for( i = 0; i < nrhs; i++ ) {
        ferr_save[i] = ferr[i];
    }
    for( i = 0; i < nrhs; i++ ) {
        berr_save[i] = berr[i];
    }

    /* Call the LAPACK routine */
    ctbrfs_( &uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx,
             ferr, berr, work, rwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_i[i] = x[i];
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
    info_i = LAPACKE_ctbrfs_work( LAPACK_COL_MAJOR, uplo_i, trans_i, diag_i,
                                  n_i, kd_i, nrhs_i, ab_i, ldab_i, b_i, ldb_i,
                                  x_i, ldx_i, ferr_i, berr_i, work_i, rwork_i );

    failed = compare_ctbrfs( ferr, ferr_i, berr, berr_i, info, info_i, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ctbrfs\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ctbrfs\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_i[i] = x[i];
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
    info_i = LAPACKE_ctbrfs( LAPACK_COL_MAJOR, uplo_i, trans_i, diag_i, n_i,
                             kd_i, nrhs_i, ab_i, ldab_i, b_i, ldb_i, x_i, ldx_i,
                             ferr_i, berr_i );

    failed = compare_ctbrfs( ferr, ferr_i, berr, berr_i, info, info_i, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ctbrfs\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ctbrfs\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_i[i] = x[i];
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

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, kd+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_ctbrfs_work( LAPACK_ROW_MAJOR, uplo_i, trans_i, diag_i,
                                  n_i, kd_i, nrhs_i, ab_r, ldab_r, b_r, ldb_r,
                                  x_r, ldx_r, ferr_i, berr_i, work_i, rwork_i );

    failed = compare_ctbrfs( ferr, ferr_i, berr, berr_i, info, info_i, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ctbrfs\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ctbrfs\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < ldb*nrhs; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldx*nrhs; i++ ) {
        x_i[i] = x[i];
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
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, kd+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_ctbrfs( LAPACK_ROW_MAJOR, uplo_i, trans_i, diag_i, n_i,
                             kd_i, nrhs_i, ab_r, ldab_r, b_r, ldb_r, x_r, ldx_r,
                             ferr_i, berr_i );

    failed = compare_ctbrfs( ferr, ferr_i, berr, berr_i, info, info_i, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ctbrfs\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ctbrfs\n" );
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

/* Auxiliary function: ctbrfs scalar parameters initialization */
static void init_scalars_ctbrfs( char *uplo, char *trans, char *diag,
                                 lapack_int *n, lapack_int *kd,
                                 lapack_int *nrhs, lapack_int *ldab,
                                 lapack_int *ldb, lapack_int *ldx )
{
    *uplo = 'L';
    *trans = 'N';
    *diag = 'N';
    *n = 4;
    *kd = 2;
    *nrhs = 2;
    *ldab = 9;
    *ldb = 8;
    *ldx = 8;

    return;
}

/* Auxiliary functions: ctbrfs array parameters initialization */
static void init_ab( lapack_int size, lapack_complex_float *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    ab[0] = lapack_make_complex_float( -1.940000057e+000, 4.429999828e+000 );
    ab[9] = lapack_make_complex_float( 4.119999886e+000, -4.269999981e+000 );
    ab[18] = lapack_make_complex_float( 4.300000072e-001, -2.660000086e+000 );
    ab[27] = lapack_make_complex_float( 4.399999976e-001, 1.000000015e-001 );
    ab[1] = lapack_make_complex_float( -3.390000105e+000, 3.440000057e+000 );
    ab[10] = lapack_make_complex_float( -1.840000033e+000, 5.530000210e+000 );
    ab[19] = lapack_make_complex_float( 1.740000010e+000, -3.999999911e-002 );
    ab[28] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[2] = lapack_make_complex_float( 1.620000005e+000, 3.680000067e+000 );
    ab[11] = lapack_make_complex_float( -2.769999981e+000, -1.929999948e+000 );
    ab[20] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[29] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_b( lapack_int size, lapack_complex_float *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    b[0] = lapack_make_complex_float( -8.859999657e+000, -3.880000114e+000 );
    b[8] = lapack_make_complex_float( -2.409000015e+001, -5.269999981e+000 );
    b[1] = lapack_make_complex_float( -1.556999969e+001, -2.340999985e+001 );
    b[9] = lapack_make_complex_float( -5.797000122e+001, 8.140000343e+000 );
    b[2] = lapack_make_complex_float( -7.630000114e+000, 2.278000069e+001 );
    b[10] = lapack_make_complex_float( 1.909000015e+001, -2.951000023e+001 );
    b[3] = lapack_make_complex_float( -1.473999977e+001, -2.400000095e+000 );
    b[11] = lapack_make_complex_float( 1.917000008e+001, 2.132999992e+001 );
}
static void init_x( lapack_int size, lapack_complex_float *x ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        x[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    x[0] = lapack_make_complex_float( 0.000000000e+000, 2.000000000e+000 );
    x[8] = lapack_make_complex_float( 1.000000119e+000, 5.000000000e+000 );
    x[1] = lapack_make_complex_float( 1.000000000e+000, -3.000000000e+000 );
    x[9] = lapack_make_complex_float( -7.000000000e+000, -1.999999642e+000 );
    x[2] = lapack_make_complex_float( -4.000000000e+000, -5.000000000e+000 );
    x[10] = lapack_make_complex_float( 2.999999285e+000, 4.000000000e+000 );
    x[3] = lapack_make_complex_float( 1.999999166e+000, -9.999992847e-001 );
    x[11] = lapack_make_complex_float( -5.999998093e+000, -8.999999046e+000 );
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

/* Auxiliary function: C interface to ctbrfs results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ctbrfs( float *ferr, float *ferr_i, float *berr,
                           float *berr_i, lapack_int info, lapack_int info_i,
                           lapack_int nrhs )
{
    lapack_int i;
    int failed = 0;
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
