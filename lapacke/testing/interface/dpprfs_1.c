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
* dpprfs_1 is the test program for the C interface to LAPACK
* routine dpprfs
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

static void init_scalars_dpprfs( char *uplo, lapack_int *n, lapack_int *nrhs,
                                 lapack_int *ldb, lapack_int *ldx );
static void init_ap( lapack_int size, double *ap );
static void init_afp( lapack_int size, double *afp );
static void init_b( lapack_int size, double *b );
static void init_x( lapack_int size, double *x );
static void init_ferr( lapack_int size, double *ferr );
static void init_berr( lapack_int size, double *berr );
static void init_work( lapack_int size, double *work );
static void init_iwork( lapack_int size, lapack_int *iwork );
static int compare_dpprfs( double *x, double *x_i, double *ferr, double *ferr_i,
                           double *berr, double *berr_i, lapack_int info,
                           lapack_int info_i, lapack_int ldx, lapack_int nrhs );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
    lapack_int n, n_i;
    lapack_int nrhs, nrhs_i;
    lapack_int ldb, ldb_i;
    lapack_int ldb_r;
    lapack_int ldx, ldx_i;
    lapack_int ldx_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *ap = NULL, *ap_i = NULL;
    double *afp = NULL, *afp_i = NULL;
    double *b = NULL, *b_i = NULL;
    double *x = NULL, *x_i = NULL;
    double *ferr = NULL, *ferr_i = NULL;
    double *berr = NULL, *berr_i = NULL;
    double *work = NULL, *work_i = NULL;
    lapack_int *iwork = NULL, *iwork_i = NULL;
    double *x_save = NULL;
    double *ferr_save = NULL;
    double *berr_save = NULL;
    double *ap_r = NULL;
    double *afp_r = NULL;
    double *b_r = NULL;
    double *x_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dpprfs( &uplo, &n, &nrhs, &ldb, &ldx );
    ldb_r = nrhs+2;
    ldx_r = nrhs+2;
    uplo_i = uplo;
    n_i = n;
    nrhs_i = nrhs;
    ldb_i = ldb;
    ldx_i = ldx;

    /* Allocate memory for the LAPACK routine arrays */
    ap = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );
    afp = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );
    b = (double *)LAPACKE_malloc( ldb*nrhs * sizeof(double) );
    x = (double *)LAPACKE_malloc( ldx*nrhs * sizeof(double) );
    ferr = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    berr = (double *)LAPACKE_malloc( nrhs * sizeof(double) );
    work = (double *)LAPACKE_malloc( 3*n * sizeof(double) );
    iwork = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    ap_i = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );
    afp_i = (double *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(double) );
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
    ap_r = (double *)LAPACKE_malloc( n*(n+1)/2 * sizeof(double) );
    afp_r = (double *)LAPACKE_malloc( n*(n+1)/2 * sizeof(double) );
    b_r = (double *)LAPACKE_malloc( n*(nrhs+2) * sizeof(double) );
    x_r = (double *)LAPACKE_malloc( n*(nrhs+2) * sizeof(double) );

    /* Initialize input arrays */
    init_ap( (n*(n+1)/2), ap );
    init_afp( (n*(n+1)/2), afp );
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
    dpprfs_( &uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work,
             iwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        afp_i[i] = afp[i];
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
    info_i = LAPACKE_dpprfs_work( LAPACK_COL_MAJOR, uplo_i, n_i, nrhs_i, ap_i,
                                  afp_i, b_i, ldb_i, x_i, ldx_i, ferr_i, berr_i,
                                  work_i, iwork_i );

    failed = compare_dpprfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dpprfs\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dpprfs\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        afp_i[i] = afp[i];
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
    info_i = LAPACKE_dpprfs( LAPACK_COL_MAJOR, uplo_i, n_i, nrhs_i, ap_i, afp_i,
                             b_i, ldb_i, x_i, ldx_i, ferr_i, berr_i );

    failed = compare_dpprfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dpprfs\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dpprfs\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        afp_i[i] = afp[i];
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

    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, n, afp_i, afp_r );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_dpprfs_work( LAPACK_ROW_MAJOR, uplo_i, n_i, nrhs_i, ap_r,
                                  afp_r, b_r, ldb_r, x_r, ldx_r, ferr_i, berr_i,
                                  work_i, iwork_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_dpprfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dpprfs\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dpprfs\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        afp_i[i] = afp[i];
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
    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    LAPACKE_dpp_trans( LAPACK_COL_MAJOR, uplo, n, afp_i, afp_r );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, nrhs, b_i, ldb, b_r, nrhs+2 );
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, nrhs, x_i, ldx, x_r, nrhs+2 );
    info_i = LAPACKE_dpprfs( LAPACK_ROW_MAJOR, uplo_i, n_i, nrhs_i, ap_r, afp_r,
                             b_r, ldb_r, x_r, ldx_r, ferr_i, berr_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, nrhs, x_r, nrhs+2, x_i, ldx );

    failed = compare_dpprfs( x, x_i, ferr, ferr_i, berr, berr_i, info, info_i,
                             ldx, nrhs );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dpprfs\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dpprfs\n" );
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
    if( afp != NULL ) {
        LAPACKE_free( afp );
    }
    if( afp_i != NULL ) {
        LAPACKE_free( afp_i );
    }
    if( afp_r != NULL ) {
        LAPACKE_free( afp_r );
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

/* Auxiliary function: dpprfs scalar parameters initialization */
static void init_scalars_dpprfs( char *uplo, lapack_int *n, lapack_int *nrhs,
                                 lapack_int *ldb, lapack_int *ldx )
{
    *uplo = 'L';
    *n = 4;
    *nrhs = 2;
    *ldb = 8;
    *ldx = 8;

    return;
}

/* Auxiliary functions: dpprfs array parameters initialization */
static void init_ap( lapack_int size, double *ap ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ap[i] = 0;
    }
    ap[0] = 4.16000000000000010e+000;
    ap[1] = -3.12000000000000010e+000;
    ap[2] = 5.60000000000000050e-001;
    ap[3] = -1.00000000000000010e-001;
    ap[4] = 5.03000000000000020e+000;
    ap[5] = -8.29999999999999960e-001;
    ap[6] = 1.17999999999999990e+000;
    ap[7] = 7.60000000000000010e-001;
    ap[8] = 3.40000000000000020e-001;
    ap[9] = 1.17999999999999990e+000;
}
static void init_afp( lapack_int size, double *afp ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        afp[i] = 0;
    }
    afp[0] = 2.03960780543711410e+000;
    afp[1] = -1.52970585407783540e+000;
    afp[2] = 2.74562589193457660e-001;
    afp[3] = -4.90290337845460110e-002;
    afp[4] = 1.64012194668567270e+000;
    afp[5] = -2.49981411948373810e-001;
    afp[6] = 6.73730390738910060e-001;
    afp[7] = 7.88748805574805310e-001;
    afp[8] = 6.61657563374256520e-001;
    afp[9] = 5.34689426929868320e-001;
}
static void init_b( lapack_int size, double *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = 0;
    }
    b[0] = 8.69999999999999930e+000;  /* b[0,0] */
    b[8] = 8.30000000000000070e+000;  /* b[0,1] */
    b[1] = -1.33500000000000000e+001;  /* b[1,0] */
    b[9] = 2.12999999999999990e+000;  /* b[1,1] */
    b[2] = 1.88999999999999990e+000;  /* b[2,0] */
    b[10] = 1.61000000000000010e+000;  /* b[2,1] */
    b[3] = -4.13999999999999970e+000;  /* b[3,0] */
    b[11] = 5.00000000000000000e+000;  /* b[3,1] */
}
static void init_x( lapack_int size, double *x ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        x[i] = 0;
    }
    x[0] = 1.00000000000000020e+000;  /* x[0,0] */
    x[8] = 3.99999999999999820e+000;  /* x[0,1] */
    x[1] = -9.99999999999998560e-001;  /* x[1,0] */
    x[9] = 2.99999999999999820e+000;  /* x[1,1] */
    x[2] = 2.00000000000000310e+000;  /* x[2,0] */
    x[10] = 1.99999999999999710e+000;  /* x[2,1] */
    x[3] = -3.00000000000000270e+000;  /* x[3,0] */
    x[11] = 1.00000000000000270e+000;  /* x[3,1] */
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

/* Auxiliary function: C interface to dpprfs results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dpprfs( double *x, double *x_i, double *ferr, double *ferr_i,
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
