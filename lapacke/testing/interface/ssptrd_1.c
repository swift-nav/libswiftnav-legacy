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
* ssptrd_1 is the test program for the C interface to LAPACK
* routine ssptrd
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

static void init_scalars_ssptrd( char *uplo, lapack_int *n );
static void init_ap( lapack_int size, float *ap );
static void init_d( lapack_int size, float *d );
static void init_e( lapack_int size, float *e );
static void init_tau( lapack_int size, float *tau );
static int compare_ssptrd( float *ap, float *ap_i, float *d, float *d_i,
                           float *e, float *e_i, float *tau, float *tau_i,
                           lapack_int info, lapack_int info_i, lapack_int n );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
    lapack_int n, n_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    float *ap = NULL, *ap_i = NULL;
    float *d = NULL, *d_i = NULL;
    float *e = NULL, *e_i = NULL;
    float *tau = NULL, *tau_i = NULL;
    float *ap_save = NULL;
    float *d_save = NULL;
    float *e_save = NULL;
    float *tau_save = NULL;
    float *ap_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ssptrd( &uplo, &n );
    uplo_i = uplo;
    n_i = n;

    /* Allocate memory for the LAPACK routine arrays */
    ap = (float *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(float) );
    d = (float *)LAPACKE_malloc( n * sizeof(float) );
    e = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );
    tau = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );

    /* Allocate memory for the C interface function arrays */
    ap_i = (float *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(float) );
    d_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    e_i = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );
    tau_i = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );

    /* Allocate memory for the backup arrays */
    ap_save = (float *)LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(float) );
    d_save = (float *)LAPACKE_malloc( n * sizeof(float) );
    e_save = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );
    tau_save = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    ap_r = (float *)LAPACKE_malloc( n*(n+1)/2 * sizeof(float) );

    /* Initialize input arrays */
    init_ap( (n*(n+1)/2), ap );
    init_d( n, d );
    init_e( (n-1), e );
    init_tau( (n-1), tau );

    /* Backup the ouptut arrays */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_save[i] = ap[i];
    }
    for( i = 0; i < n; i++ ) {
        d_save[i] = d[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_save[i] = e[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        tau_save[i] = tau[i];
    }

    /* Call the LAPACK routine */
    ssptrd_( &uplo, &n, ap, d, e, tau, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        tau_i[i] = tau_save[i];
    }
    info_i = LAPACKE_ssptrd_work( LAPACK_COL_MAJOR, uplo_i, n_i, ap_i, d_i, e_i,
                                  tau_i );

    failed = compare_ssptrd( ap, ap_i, d, d_i, e, e_i, tau, tau_i, info, info_i,
                             n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ssptrd\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ssptrd\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        tau_i[i] = tau_save[i];
    }
    info_i = LAPACKE_ssptrd( LAPACK_COL_MAJOR, uplo_i, n_i, ap_i, d_i, e_i,
                             tau_i );

    failed = compare_ssptrd( ap, ap_i, d, d_i, e, e_i, tau, tau_i, info, info_i,
                             n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ssptrd\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ssptrd\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        tau_i[i] = tau_save[i];
    }

    LAPACKE_spp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    info_i = LAPACKE_ssptrd_work( LAPACK_ROW_MAJOR, uplo_i, n_i, ap_r, d_i, e_i,
                                  tau_i );

    LAPACKE_spp_trans( LAPACK_ROW_MAJOR, uplo, n, ap_r, ap_i );

    failed = compare_ssptrd( ap, ap_i, d, d_i, e, e_i, tau, tau_i, info, info_i,
                             n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ssptrd\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ssptrd\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        tau_i[i] = tau_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_spp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    info_i = LAPACKE_ssptrd( LAPACK_ROW_MAJOR, uplo_i, n_i, ap_r, d_i, e_i,
                             tau_i );

    LAPACKE_spp_trans( LAPACK_ROW_MAJOR, uplo, n, ap_r, ap_i );

    failed = compare_ssptrd( ap, ap_i, d, d_i, e, e_i, tau, tau_i, info, info_i,
                             n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ssptrd\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ssptrd\n" );
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
    if( ap_save != NULL ) {
        LAPACKE_free( ap_save );
    }
    if( d != NULL ) {
        LAPACKE_free( d );
    }
    if( d_i != NULL ) {
        LAPACKE_free( d_i );
    }
    if( d_save != NULL ) {
        LAPACKE_free( d_save );
    }
    if( e != NULL ) {
        LAPACKE_free( e );
    }
    if( e_i != NULL ) {
        LAPACKE_free( e_i );
    }
    if( e_save != NULL ) {
        LAPACKE_free( e_save );
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

    return 0;
}

/* Auxiliary function: ssptrd scalar parameters initialization */
static void init_scalars_ssptrd( char *uplo, lapack_int *n )
{
    *uplo = 'L';
    *n = 4;

    return;
}

/* Auxiliary functions: ssptrd array parameters initialization */
static void init_ap( lapack_int size, float *ap ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ap[i] = 0;
    }
    ap[0] = 2.069999933e+000;
    ap[1] = 3.869999886e+000;
    ap[2] = 4.199999809e+000;
    ap[3] = -1.149999976e+000;
    ap[4] = -2.099999934e-001;
    ap[5] = 1.870000005e+000;
    ap[6] = 6.299999952e-001;
    ap[7] = 1.149999976e+000;
    ap[8] = 2.059999943e+000;
    ap[9] = -1.809999943e+000;
}
static void init_d( lapack_int size, float *d ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        d[i] = 0;
    }
}
static void init_e( lapack_int size, float *e ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        e[i] = 0;
    }
}
static void init_tau( lapack_int size, float *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = 0;
    }
}

/* Auxiliary function: C interface to ssptrd results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ssptrd( float *ap, float *ap_i, float *d, float *d_i,
                           float *e, float *e_i, float *tau, float *tau_i,
                           lapack_int info, lapack_int info_i, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        failed += compare_floats(ap[i],ap_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_floats(d[i],d_i[i]);
    }
    for( i = 0; i < (n-1); i++ ) {
        failed += compare_floats(e[i],e_i[i]);
    }
    for( i = 0; i < (n-1); i++ ) {
        failed += compare_floats(tau[i],tau_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
