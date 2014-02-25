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
* sgbrfs_1 is the test program for the C interface to LAPACK
* routine sgbrfs
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

static void init_scalars_sgbrfs( char *trans, lapack_int *n, lapack_int *kl,
                                 lapack_int *ku, lapack_int *nrhs,
                                 lapack_int *ldab, lapack_int *ldafb,
                                 lapack_int *ldb, lapack_int *ldx );
static void init_ab( lapack_int size, float *ab );
static void init_afb( lapack_int size, float *afb );
static void init_ipiv( lapack_int size, lapack_int *ipiv );
static void init_b( lapack_int size, float *b );
static void init_x( lapack_int size, float *x );
static void init_ferr( lapack_int size, float *ferr );
static void init_berr( lapack_int size, float *berr );
static void init_work( lapack_int size, float *work );
static void init_iwork( lapack_int size, lapack_int *iwork );
static int compare_sgbrfs( float *x, float *x_i, float *ferr, float *ferr_i,
                           float *berr, float *berr_i, lapack_int info,
                           lapack_int info_i, lapack_int ldx, lapack_int nrhs );

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
    float *ab = NULL, *ab_i = NULL;
    float *afb = NULL, *afb_i = NULL;
    lapack_int *ipiv = NULL, *ipiv_i = NULL;
    float *b = NULL, *b_i = NULL;
    float *x = NULL, *x_i = NULL;
    float *ferr = NULL, *ferr_i = NULL;
    float *berr = NULL, *berr_i = NULL;
    float *work = NULL, *work_i = NULL;
    lapack_int *iwork = NULL, *iwork_i = NULL;
    float *x_save = NULL;
    float *ferr_save = NULL;
    float *berr_save = NULL;
    float *ab_r = NULL;
    float *afb_r = NULL;
    float *b_r = NULL;
    float *x_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_sgbrfs( &trans, &n, &kl, &ku, &nrhs, &ldab, &ldafb, &ldb,
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
    ab = (float *)LAPACKE_malloc( ldab*n * sizeof(float) );
    afb = (float *)LAPACKE_malloc( ldafb*n * sizeof(float) );
    ipiv = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b = (float *)LAPACKE_malloc( ldb*nrhs * sizeof(float) );
    x = (float *)LAPACKE_malloc( ldx*nrhs * sizeof(float) );
    ferr = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    berr = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    work = (float *)LAPACKE_malloc( 3*n * sizeof(float) );
    iwork = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (float *)LAPACKE_malloc( ldab*n * sizeof(float) );
    afb_i = (float *)LAPACKE_malloc( ldafb*n * sizeof(float) );
    ipiv_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    b_i = (float *)LAPACKE_malloc( ldb*nrhs * sizeof(float) );
    x_i = (float *)LAPACKE_malloc( ldx*nrhs * sizeof(float) );
    ferr_i = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    berr_i = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    work_i = (float *)LAPACKE_malloc( 3*n * sizeof(float) );
    iwork_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    x_save = (float *)LAPACKE_malloc( ldx*nrhs * sizeof(float) );
    ferr_save = (float *)LAPACKE_malloc( nrhs * sizeof(float) );
    berr_save = (float *)LAPACKE_malloc( nrhs * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    ab_r = (float *)LAPACKE_malloc( (kl+ku+1)*(n+2) * sizeof(float) );
    afb_r = (float *)LAPACKE_malloc( ((2*kl+ku+1)*(n+2)) * sizeof(float) );
    b_r = (float *)LAPACKE_malloc( n*(nrhs+2) * sizeof(float) );
    x_r = (float *)LAPACKE_malloc( n*(nrhs+2) * sizeof(float) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_afb( ldafb*n, afb );
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
    sgbrfs_( &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb,
             x, &ldx, ferr, berr, work, iwork, &info );

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
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_sgbrfs_work( LAPACK_COL_MAJOR, trans_i, n_i, kl_i, ku_i,
                                  nrhs_i, ab_i, ldab_i, afb_i, ldafb_i, ipiv_i,
                                  b_i, ldb_i, x_i, ldx_i, ferr_i, berr_i,
                                  work_i, iwork_i );

    failed = compare_sgbrfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to sgbrfs\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to sgbrfs\n" );
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
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_sgbrfs( LAPACK_COL_MAJOR, trans_i, n_i, kl_i, ku_i, nrhs_i,
                             ab_i, ldab_i, afb_i, ldafb_i, ipiv_i, b_i, ldb_i,
                             x_i, ldx_i, ferr_i, berr_i );

    failed = compare_sgbrfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to sgbrfs\n" );
    } else {
        printf( "FAILED: column-major high-level interface to sgbrfs\n" );
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
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }

    LAPACKE_sge_trans( LAPACK_COL_MAJOR, kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, 2*kl+ku+1, n, afb_i, ldafb, afb_r,
                       n+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_sgbrfs_work( LAPACK_ROW_MAJOR, trans_i, n_i, kl_i, ku_i,
                                  nrhs_i, ab_r, ldab_r, afb_r, ldafb_r, ipiv_i,
                                  b_r, ldb_r, x_r, ldx_r, ferr_i, berr_i,
                                  work_i, iwork_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_sgbrfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to sgbrfs\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to sgbrfs\n" );
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
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, 2*kl+ku+1, n, afb_i, ldafb, afb_r,
                       n+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_sgbrfs( LAPACK_ROW_MAJOR, trans_i, n_i, kl_i, ku_i, nrhs_i,
                             ab_r, ldab_r, afb_r, ldafb_r, ipiv_i, b_r, ldb_r,
                             x_r, ldx_r, ferr_i, berr_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_sgbrfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to sgbrfs\n" );
    } else {
        printf( "FAILED: row-major high-level interface to sgbrfs\n" );
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
    if( iwork != NULL ) {
        LAPACKE_free( iwork );
    }
    if( iwork_i != NULL ) {
        LAPACKE_free( iwork_i );
    }

    return 0;
}

/* Auxiliary function: sgbrfs scalar parameters initialization */
static void init_scalars_sgbrfs( char *trans, lapack_int *n, lapack_int *kl,
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

/* Auxiliary functions: sgbrfs array parameters initialization */
static void init_ab( lapack_int size, float *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = 0;
    }
    ab[0] = 0.000000000e+000;  /* ab[0,0] */
    ab[17] = 0.000000000e+000;  /* ab[0,1] */
    ab[34] = -3.660000086e+000;  /* ab[0,2] */
    ab[51] = -2.130000114e+000;  /* ab[0,3] */
    ab[1] = 0.000000000e+000;  /* ab[1,0] */
    ab[18] = 2.539999962e+000;  /* ab[1,1] */
    ab[35] = -2.730000019e+000;  /* ab[1,2] */
    ab[52] = 4.070000172e+000;  /* ab[1,3] */
    ab[2] = -2.300000042e-001;  /* ab[2,0] */
    ab[19] = 2.460000038e+000;  /* ab[2,1] */
    ab[36] = 2.460000038e+000;  /* ab[2,2] */
    ab[53] = -3.819999933e+000;  /* ab[2,3] */
    ab[3] = -6.980000019e+000;  /* ab[3,0] */
    ab[20] = 2.559999943e+000;  /* ab[3,1] */
    ab[37] = -4.780000210e+000;  /* ab[3,2] */
    ab[54] = 0.000000000e+000;  /* ab[3,3] */
}
static void init_afb( lapack_int size, float *afb ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        afb[i] = 0;
    }
    afb[0] = 0.000000000e+000;  /* afb[0,0] */
    afb[25] = 0.000000000e+000;  /* afb[0,1] */
    afb[50] = 0.000000000e+000;  /* afb[0,2] */
    afb[75] = -2.130000114e+000;  /* afb[0,3] */
    afb[1] = 0.000000000e+000;  /* afb[1,0] */
    afb[26] = 0.000000000e+000;  /* afb[1,1] */
    afb[51] = -2.730000019e+000;  /* afb[1,2] */
    afb[76] = 4.070000172e+000;  /* afb[1,3] */
    afb[2] = 0.000000000e+000;  /* afb[2,0] */
    afb[27] = 2.460000038e+000;  /* afb[2,1] */
    afb[52] = 2.460000038e+000;  /* afb[2,2] */
    afb[77] = -3.839143991e+000;  /* afb[2,3] */
    afb[3] = -6.980000019e+000;  /* afb[3,0] */
    afb[28] = 2.559999943e+000;  /* afb[3,1] */
    afb[53] = -5.932930470e+000;  /* afb[3,2] */
    afb[78] = -7.269061804e-001;  /* afb[3,3] */
    afb[4] = 3.295128793e-002;  /* afb[4,0] */
    afb[29] = 9.605233669e-001;  /* afb[4,1] */
    afb[54] = 8.056727648e-001;  /* afb[4,2] */
    afb[79] = 0.000000000e+000;  /* afb[4,3] */
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
static void init_b( lapack_int size, float *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = 0;
    }
    b[0] = 4.420000076e+000;  /* b[0,0] */
    b[8] = -3.600999832e+001;  /* b[0,1] */
    b[1] = 2.712999916e+001;  /* b[1,0] */
    b[9] = -3.167000008e+001;  /* b[1,1] */
    b[2] = -6.139999866e+000;  /* b[2,0] */
    b[10] = -1.159999967e+000;  /* b[2,1] */
    b[3] = 1.050000000e+001;  /* b[3,0] */
    b[11] = -2.581999969e+001;  /* b[3,1] */
}
static void init_x( lapack_int size, float *x ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        x[i] = 0;
    }
    x[0] = -1.999999285e+000;  /* x[0,0] */
    x[8] = 1.000001073e+000;  /* x[0,1] */
    x[1] = 3.000001669e+000;  /* x[1,0] */
    x[9] = -3.999997854e+000;  /* x[1,1] */
    x[2] = 1.000001073e+000;  /* x[2,0] */
    x[10] = 7.000000954e+000;  /* x[2,1] */
    x[3] = -4.000001431e+000;  /* x[3,0] */
    x[11] = -2.000002146e+000;  /* x[3,1] */
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

/* Auxiliary function: C interface to sgbrfs results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_sgbrfs( float *x, float *x_i, float *ferr, float *ferr_i,
                           float *berr, float *berr_i, lapack_int info,
                           lapack_int info_i, lapack_int ldx, lapack_int nrhs )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldx*nrhs; i++ ) {
        failed += compare_floats(x[i],x_i[i]);
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
