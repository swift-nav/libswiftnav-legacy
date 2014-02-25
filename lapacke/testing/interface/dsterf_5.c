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
* dsterf_5 is the test program for the C interface to LAPACK
* routine dsterf
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

static void init_scalars_dsterf( lapack_int *n );
static void init_d( lapack_int size, double *d );
static void init_e( lapack_int size, double *e );
static int compare_dsterf( double *d, double *d_i, double *e, double *e_i,
                           lapack_int info, lapack_int info_i, lapack_int n );

int main(void)
{
    /* Local scalars */
    lapack_int n, n_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *d = NULL, *d_i = NULL;
    double *e = NULL, *e_i = NULL;
    double *d_save = NULL;
    double *e_save = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dsterf( &n );
    n_i = n;

    /* Allocate memory for the LAPACK routine arrays */
    d = (double *)LAPACKE_malloc( n * sizeof(double) );
    e = (double *)LAPACKE_malloc( (n-1) * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    d_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    e_i = (double *)LAPACKE_malloc( (n-1) * sizeof(double) );

    /* Allocate memory for the backup arrays */
    d_save = (double *)LAPACKE_malloc( n * sizeof(double) );
    e_save = (double *)LAPACKE_malloc( (n-1) * sizeof(double) );

    /* Allocate memory for the row-major arrays */

    /* Initialize input arrays */
    init_d( n, d );
    init_e( (n-1), e );

    /* Backup the ouptut arrays */
    for( i = 0; i < n; i++ ) {
        d_save[i] = d[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_save[i] = e[i];
    }

    /* Call the LAPACK routine */
    dsterf_( &n, d, e, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    info_i = LAPACKE_dsterf_work( n_i, d_i, e_i );

    failed = compare_dsterf( d, d_i, e, e_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dsterf\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dsterf\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    info_i = LAPACKE_dsterf( n_i, d_i, e_i );

    failed = compare_dsterf( d, d_i, e, e_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dsterf\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dsterf\n" );
    }

    failed = compare_dsterf( d, d_i, e, e_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dsterf\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dsterf\n" );
    }

    failed = compare_dsterf( d, d_i, e, e_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dsterf\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dsterf\n" );
    }

    /* Release memory */
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

    return 0;
}

/* Auxiliary function: dsterf scalar parameters initialization */
static void init_scalars_dsterf( lapack_int *n )
{
    *n = 4;

    return;
}

/* Auxiliary functions: dsterf array parameters initialization */
static void init_d( lapack_int size, double *d ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        d[i] = 0;
    }
    d[0] = -2.27863777089783250e+000;
    d[1] = -1.33501501572863330e-001;
    d[2] = -1.46357763012886500e-001;
    d[3] = -1.93038324810632100e+000;
}
static void init_e( lapack_int size, double *e ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        e[i] = 0;
    }
    e[0] = -4.33774701789705920e+000;
    e[1] = -2.02251863394581430e+000;
    e[2] = -1.79232441360670490e+000;
}

/* Auxiliary function: C interface to dsterf results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dsterf( double *d, double *d_i, double *e, double *e_i,
                           lapack_int info, lapack_int info_i, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < n; i++ ) {
        failed += compare_doubles(d[i],d_i[i]);
    }
    for( i = 0; i < (n-1); i++ ) {
        failed += compare_doubles(e[i],e_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
