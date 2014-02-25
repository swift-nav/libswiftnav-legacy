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
* cbdsqr_2 is the test program for the C interface to LAPACK
* routine cbdsqr
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

static void init_scalars_cbdsqr( char *uplo, lapack_int *n, lapack_int *ncvt,
                                 lapack_int *nru, lapack_int *ncc,
                                 lapack_int *ldvt, lapack_int *ldu,
                                 lapack_int *ldc );
static void init_d( lapack_int size, float *d );
static void init_e( lapack_int size, float *e );
static void init_vt( lapack_int size, lapack_complex_float *vt );
static void init_u( lapack_int size, lapack_complex_float *u );
static void init_c( lapack_int size, lapack_complex_float *c );
static void init_work( lapack_int size, float *work );
static int compare_cbdsqr( float *d, float *d_i, float *e, float *e_i,
                           lapack_complex_float *vt, lapack_complex_float *vt_i,
                           lapack_complex_float *u, lapack_complex_float *u_i,
                           lapack_complex_float *c, lapack_complex_float *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int ldu, lapack_int ldvt, lapack_int n,
                           lapack_int ncc, lapack_int ncvt, lapack_int nru );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
    lapack_int n, n_i;
    lapack_int ncvt, ncvt_i;
    lapack_int nru, nru_i;
    lapack_int ncc, ncc_i;
    lapack_int ldvt, ldvt_i;
    lapack_int ldvt_r;
    lapack_int ldu, ldu_i;
    lapack_int ldu_r;
    lapack_int ldc, ldc_i;
    lapack_int ldc_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    float *d = NULL, *d_i = NULL;
    float *e = NULL, *e_i = NULL;
    lapack_complex_float *vt = NULL, *vt_i = NULL;
    lapack_complex_float *u = NULL, *u_i = NULL;
    lapack_complex_float *c = NULL, *c_i = NULL;
    float *work = NULL, *work_i = NULL;
    float *d_save = NULL;
    float *e_save = NULL;
    lapack_complex_float *vt_save = NULL;
    lapack_complex_float *u_save = NULL;
    lapack_complex_float *c_save = NULL;
    lapack_complex_float *vt_r = NULL;
    lapack_complex_float *u_r = NULL;
    lapack_complex_float *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_cbdsqr( &uplo, &n, &ncvt, &nru, &ncc, &ldvt, &ldu, &ldc );
    ldvt_r = ncvt+2;
    ldu_r = n+2;
    ldc_r = ncc+2;
    uplo_i = uplo;
    n_i = n;
    ncvt_i = ncvt;
    nru_i = nru;
    ncc_i = ncc;
    ldvt_i = ldvt;
    ldu_i = ldu;
    ldc_i = ldc;

    /* Allocate memory for the LAPACK routine arrays */
    d = (float *)LAPACKE_malloc( n * sizeof(float) );
    e = (float *)LAPACKE_malloc( n * sizeof(float) );
    vt = (lapack_complex_float *)
        LAPACKE_malloc( ldvt*ncvt * sizeof(lapack_complex_float) );
    u = (lapack_complex_float *)
        LAPACKE_malloc( ldu*n * sizeof(lapack_complex_float) );
    c = (lapack_complex_float *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_float) );
    work = (float *)LAPACKE_malloc( 4*n * sizeof(float) );

    /* Allocate memory for the C interface function arrays */
    d_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    e_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    vt_i = (lapack_complex_float *)
        LAPACKE_malloc( ldvt*ncvt * sizeof(lapack_complex_float) );
    u_i = (lapack_complex_float *)
        LAPACKE_malloc( ldu*n * sizeof(lapack_complex_float) );
    c_i = (lapack_complex_float *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_float) );
    work_i = (float *)LAPACKE_malloc( 4*n * sizeof(float) );

    /* Allocate memory for the backup arrays */
    d_save = (float *)LAPACKE_malloc( n * sizeof(float) );
    e_save = (float *)LAPACKE_malloc( n * sizeof(float) );
    vt_save = (lapack_complex_float *)
        LAPACKE_malloc( ldvt*ncvt * sizeof(lapack_complex_float) );
    u_save = (lapack_complex_float *)
        LAPACKE_malloc( ldu*n * sizeof(lapack_complex_float) );
    c_save = (lapack_complex_float *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    vt_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(ncvt+2) * sizeof(lapack_complex_float) );
    u_r = (lapack_complex_float *)
        LAPACKE_malloc( nru*(n+2) * sizeof(lapack_complex_float) );
    c_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(ncc+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_d( n, d );
    init_e( n, e );
    init_vt( ldvt*ncvt, vt );
    init_u( ldu*n, u );
    init_c( ldc*ncc, c );
    init_work( 4*n, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < n; i++ ) {
        d_save[i] = d[i];
    }
    for( i = 0; i < n; i++ ) {
        e_save[i] = e[i];
    }
    for( i = 0; i < ldvt*ncvt; i++ ) {
        vt_save[i] = vt[i];
    }
    for( i = 0; i < ldu*n; i++ ) {
        u_save[i] = u[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    cbdsqr_( &uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc,
             work, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < n; i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldvt*ncvt; i++ ) {
        vt_i[i] = vt_save[i];
    }
    for( i = 0; i < ldu*n; i++ ) {
        u_i[i] = u_save[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < 4*n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_cbdsqr_work( LAPACK_COL_MAJOR, uplo_i, n_i, ncvt_i, nru_i,
                                  ncc_i, d_i, e_i, vt_i, ldvt_i, u_i, ldu_i,
                                  c_i, ldc_i, work_i );

    failed = compare_cbdsqr( d, d_i, e, e_i, vt, vt_i, u, u_i, c, c_i, info,
                             info_i, ldc, ldu, ldvt, n, ncc, ncvt, nru );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to cbdsqr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to cbdsqr\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < n; i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldvt*ncvt; i++ ) {
        vt_i[i] = vt_save[i];
    }
    for( i = 0; i < ldu*n; i++ ) {
        u_i[i] = u_save[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < 4*n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_cbdsqr( LAPACK_COL_MAJOR, uplo_i, n_i, ncvt_i, nru_i,
                             ncc_i, d_i, e_i, vt_i, ldvt_i, u_i, ldu_i, c_i,
                             ldc_i );

    failed = compare_cbdsqr( d, d_i, e, e_i, vt, vt_i, u, u_i, c, c_i, info,
                             info_i, ldc, ldu, ldvt, n, ncc, ncvt, nru );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to cbdsqr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to cbdsqr\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < n; i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldvt*ncvt; i++ ) {
        vt_i[i] = vt_save[i];
    }
    for( i = 0; i < ldu*n; i++ ) {
        u_i[i] = u_save[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < 4*n; i++ ) {
        work_i[i] = work[i];
    }

    if( ncvt != 0 ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, ncvt, vt_i, ldvt, vt_r,
                           ncvt+2 );
    }
    if( nru != 0 ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, nru, n, u_i, ldu, u_r, n+2 );
    }
    if( ncc != 0 ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, ncc, c_i, ldc, c_r, ncc+2 );
    }
    info_i = LAPACKE_cbdsqr_work( LAPACK_ROW_MAJOR, uplo_i, n_i, ncvt_i, nru_i,
                                  ncc_i, d_i, e_i, vt_r, ldvt_r, u_r, ldu_r,
                                  c_r, ldc_r, work_i );

    if( ncvt != 0 ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, ncvt, vt_r, ncvt+2, vt_i,
                           ldvt );
    }
    if( nru != 0 ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, nru, n, u_r, n+2, u_i, ldu );
    }
    if( ncc != 0 ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, ncc, c_r, ncc+2, c_i, ldc );
    }

    failed = compare_cbdsqr( d, d_i, e, e_i, vt, vt_i, u, u_i, c, c_i, info,
                             info_i, ldc, ldu, ldvt, n, ncc, ncvt, nru );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to cbdsqr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to cbdsqr\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < n; i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldvt*ncvt; i++ ) {
        vt_i[i] = vt_save[i];
    }
    for( i = 0; i < ldu*n; i++ ) {
        u_i[i] = u_save[i];
    }
    for( i = 0; i < ldc*ncc; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < 4*n; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    if( ncvt != 0 ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, ncvt, vt_i, ldvt, vt_r,
                           ncvt+2 );
    }
    if( nru != 0 ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, nru, n, u_i, ldu, u_r, n+2 );
    }
    if( ncc != 0 ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, ncc, c_i, ldc, c_r, ncc+2 );
    }
    info_i = LAPACKE_cbdsqr( LAPACK_ROW_MAJOR, uplo_i, n_i, ncvt_i, nru_i,
                             ncc_i, d_i, e_i, vt_r, ldvt_r, u_r, ldu_r, c_r,
                             ldc_r );

    if( ncvt != 0 ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, ncvt, vt_r, ncvt+2, vt_i,
                           ldvt );
    }
    if( nru != 0 ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, nru, n, u_r, n+2, u_i, ldu );
    }
    if( ncc != 0 ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, ncc, c_r, ncc+2, c_i, ldc );
    }

    failed = compare_cbdsqr( d, d_i, e, e_i, vt, vt_i, u, u_i, c, c_i, info,
                             info_i, ldc, ldu, ldvt, n, ncc, ncvt, nru );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to cbdsqr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to cbdsqr\n" );
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
    if( vt != NULL ) {
        LAPACKE_free( vt );
    }
    if( vt_i != NULL ) {
        LAPACKE_free( vt_i );
    }
    if( vt_r != NULL ) {
        LAPACKE_free( vt_r );
    }
    if( vt_save != NULL ) {
        LAPACKE_free( vt_save );
    }
    if( u != NULL ) {
        LAPACKE_free( u );
    }
    if( u_i != NULL ) {
        LAPACKE_free( u_i );
    }
    if( u_r != NULL ) {
        LAPACKE_free( u_r );
    }
    if( u_save != NULL ) {
        LAPACKE_free( u_save );
    }
    if( c != NULL ) {
        LAPACKE_free( c );
    }
    if( c_i != NULL ) {
        LAPACKE_free( c_i );
    }
    if( c_r != NULL ) {
        LAPACKE_free( c_r );
    }
    if( c_save != NULL ) {
        LAPACKE_free( c_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: cbdsqr scalar parameters initialization */
static void init_scalars_cbdsqr( char *uplo, lapack_int *n, lapack_int *ncvt,
                                 lapack_int *nru, lapack_int *ncc,
                                 lapack_int *ldvt, lapack_int *ldu,
                                 lapack_int *ldc )
{
    *uplo = 'L';
    *n = 3;
    *ncvt = 4;
    *nru = 3;
    *ncc = 0;
    *ldvt = 8;
    *ldu = 8;
    *ldc = 8;

    return;
}

/* Auxiliary functions: cbdsqr array parameters initialization */
static void init_d( lapack_int size, float *d ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        d[i] = 0;
    }
    d[0] = -2.225511312e+000;
    d[1] = 1.604558587e+000;
    d[2] = -1.673056483e+000;
}
static void init_e( lapack_int size, float *e ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        e[i] = 0;
    }
    e[0] = -1.634884596e+000;
    e[1] = 9.919312596e-001;
    e[2] = 0.000000000e+000;
}
static void init_vt( lapack_int size, lapack_complex_float *vt ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vt[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    vt[0] = lapack_make_complex_float( -1.258138418e-001, 1.617605835e-001 );
    vt[8] = lapack_make_complex_float( -2.246674895e-001, 3.864280879e-001 );
    vt[16] = lapack_make_complex_float( 3.459879458e-001, 2.156807780e-001 );
    vt[24] = lapack_make_complex_float( -7.099492550e-001, -2.965611219e-001 );
    vt[1] = lapack_make_complex_float( 4.147850871e-001, 1.795049161e-001 );
    vt[9] = lapack_make_complex_float( 1.367822587e-001, -3.976195157e-001 );
    vt[17] = lapack_make_complex_float( 6.884793043e-001, 3.386386931e-001 );
    vt[25] = lapack_make_complex_float( 1.667285562e-001, -4.942050576e-002 );
    vt[2] = lapack_make_complex_float( 4.575184286e-001, -4.806572199e-001 );
    vt[10] = lapack_make_complex_float( -2.733266354e-001, 4.980674982e-001 );
    vt[18] = lapack_make_complex_float( -2.296701074e-002, 3.860528469e-001 );
    vt[26] = lapack_make_complex_float( 1.730471253e-001, 2.394939810e-001 );
}
static void init_u( lapack_int size, lapack_complex_float *u ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        u[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    u[0] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
    u[8] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    u[16] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    u[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    u[9] = lapack_make_complex_float( -5.020265579e-001, -7.575196028e-001 );
    u[17] = lapack_make_complex_float( -3.596866727e-001, 2.115632296e-001 );
    u[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    u[10] = lapack_make_complex_float( -6.321370602e-004, 4.172925353e-001 );
    u[18] = lapack_make_complex_float( -9.073931575e-001, -5.004292727e-002 );
}
static void init_c( lapack_int size, lapack_complex_float *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}
static void init_work( lapack_int size, float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to cbdsqr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_cbdsqr( float *d, float *d_i, float *e, float *e_i,
                           lapack_complex_float *vt, lapack_complex_float *vt_i,
                           lapack_complex_float *u, lapack_complex_float *u_i,
                           lapack_complex_float *c, lapack_complex_float *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int ldu, lapack_int ldvt, lapack_int n,
                           lapack_int ncc, lapack_int ncvt, lapack_int nru )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < n; i++ ) {
        failed += compare_floats(d[i],d_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_floats(e[i],e_i[i]);
    }
    if( ncvt != 0 ) {
        for( i = 0; i < ldvt*ncvt; i++ ) {
            failed += compare_complex_floats(vt[i],vt_i[i]);
        }
    }
    if( nru != 0 ) {
        for( i = 0; i < ldu*n; i++ ) {
            failed += compare_complex_floats(u[i],u_i[i]);
        }
    }
    if( ncc != 0 ) {
        for( i = 0; i < ldc*ncc; i++ ) {
            failed += compare_complex_floats(c[i],c_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
