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
* dgbcon_1 is the test program for the C interface to LAPACK
* routine dgbcon
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

static void init_scalars_dgbcon( char *norm, lapack_int *n, lapack_int *kl,
                                 lapack_int *ku, lapack_int *ldab,
                                 double *anorm );
static void init_ab( lapack_int size, double *ab );
static void init_ipiv( lapack_int size, lapack_int *ipiv );
static void init_work( lapack_int size, double *work );
static void init_iwork( lapack_int size, lapack_int *iwork );
static int compare_dgbcon( double rcond, double rcond_i, lapack_int info,
                           lapack_int info_i );

int main(void)
{
    /* Local scalars */
    char norm, norm_i;
    lapack_int n, n_i;
    lapack_int kl, kl_i;
    lapack_int ku, ku_i;
    lapack_int ldab, ldab_i;
    lapack_int ldab_r;
    double anorm, anorm_i;
    double rcond, rcond_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *ab = NULL, *ab_i = NULL;
    lapack_int *ipiv = NULL, *ipiv_i = NULL;
    double *work = NULL, *work_i = NULL;
    lapack_int *iwork = NULL, *iwork_i = NULL;
    double *ab_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dgbcon( &norm, &n, &kl, &ku, &ldab, &anorm );
    ldab_r = n+2;
    norm_i = norm;
    n_i = n;
    kl_i = kl;
    ku_i = ku;
    ldab_i = ldab;
    anorm_i = anorm;

    /* Allocate memory for the LAPACK routine arrays */
    ab = (double *)LAPACKE_malloc( ldab*n * sizeof(double) );
    ipiv = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    work = (double *)LAPACKE_malloc( 3*n * sizeof(double) );
    iwork = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (double *)LAPACKE_malloc( ldab*n * sizeof(double) );
    ipiv_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    work_i = (double *)LAPACKE_malloc( 3*n * sizeof(double) );
    iwork_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );

    /* Allocate memory for the row-major arrays */
    ab_r = (double *)LAPACKE_malloc( ((2*kl+ku+1)*(n+2)) * sizeof(double) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_ipiv( n, ipiv );
    init_work( 3*n, work );
    init_iwork( n, iwork );

    /* Call the LAPACK routine */
    dgbcon_( &norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, &rcond, work, iwork,
             &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_dgbcon_work( LAPACK_COL_MAJOR, norm_i, n_i, kl_i, ku_i,
                                  ab_i, ldab_i, ipiv_i, anorm_i, &rcond_i,
                                  work_i, iwork_i );

    failed = compare_dgbcon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dgbcon\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dgbcon\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    info_i = LAPACKE_dgbcon( LAPACK_COL_MAJOR, norm_i, n_i, kl_i, ku_i, ab_i,
                             ldab_i, ipiv_i, anorm_i, &rcond_i );

    failed = compare_dgbcon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dgbcon\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dgbcon\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, 2*kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    info_i = LAPACKE_dgbcon_work( LAPACK_ROW_MAJOR, norm_i, n_i, kl_i, ku_i,
                                  ab_r, ldab_r, ipiv_i, anorm_i, &rcond_i,
                                  work_i, iwork_i );

    failed = compare_dgbcon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dgbcon\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dgbcon\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        ipiv_i[i] = ipiv[i];
    }
    for( i = 0; i < 3*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, 2*kl+ku+1, n, ab_i, ldab, ab_r, n+2 );
    info_i = LAPACKE_dgbcon( LAPACK_ROW_MAJOR, norm_i, n_i, kl_i, ku_i, ab_r,
                             ldab_r, ipiv_i, anorm_i, &rcond_i );

    failed = compare_dgbcon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dgbcon\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dgbcon\n" );
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

/* Auxiliary function: dgbcon scalar parameters initialization */
static void init_scalars_dgbcon( char *norm, lapack_int *n, lapack_int *kl,
                                 lapack_int *ku, lapack_int *ldab,
                                 double *anorm )
{
    *norm = '1';
    *n = 4;
    *kl = 1;
    *ku = 2;
    *ldab = 25;
    *anorm = 1.36300000000000030e+001;

    return;
}

/* Auxiliary functions: dgbcon array parameters initialization */
static void init_ab( lapack_int size, double *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = 0;
    }
    ab[0] = 0.00000000000000000e+000;  /* ab[0,0] */
    ab[25] = 0.00000000000000000e+000;  /* ab[0,1] */
    ab[50] = 0.00000000000000000e+000;  /* ab[0,2] */
    ab[75] = -2.12999999999999990e+000;  /* ab[0,3] */
    ab[1] = 0.00000000000000000e+000;  /* ab[1,0] */
    ab[26] = 0.00000000000000000e+000;  /* ab[1,1] */
    ab[51] = -2.73000000000000000e+000;  /* ab[1,2] */
    ab[76] = 4.07000000000000030e+000;  /* ab[1,3] */
    ab[2] = 0.00000000000000000e+000;  /* ab[2,0] */
    ab[27] = 2.46000000000000000e+000;  /* ab[2,1] */
    ab[52] = 2.46000000000000000e+000;  /* ab[2,2] */
    ab[77] = -3.83914387088108940e+000;  /* ab[2,3] */
    ab[3] = -6.98000000000000040e+000;  /* ab[3,0] */
    ab[28] = 2.56000000000000010e+000;  /* ab[3,1] */
    ab[53] = -5.93293047098853950e+000;  /* ab[3,2] */
    ab[78] = -7.26906663992311850e-001;  /* ab[3,3] */
    ab[4] = 3.29512893982807990e-002;  /* ab[4,0] */
    ab[29] = 9.60523370343839610e-001;  /* ab[4,1] */
    ab[54] = 8.05672681211037410e-001;  /* ab[4,2] */
    ab[79] = 0.00000000000000000e+000;  /* ab[4,3] */
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

/* Auxiliary function: C interface to dgbcon results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dgbcon( double rcond, double rcond_i, lapack_int info,
                           lapack_int info_i )
{
    int failed = 0;
    failed += compare_doubles(rcond,rcond_i);
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
