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
* zbdsqr_2 is the test program for the C interface to LAPACK
* routine zbdsqr
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

static void init_scalars_zbdsqr( char *uplo, lapack_int *n, lapack_int *ncvt,
                                 lapack_int *nru, lapack_int *ncc,
                                 lapack_int *ldvt, lapack_int *ldu,
                                 lapack_int *ldc );
static void init_d( lapack_int size, double *d );
static void init_e( lapack_int size, double *e );
static void init_vt( lapack_int size, lapack_complex_double *vt );
static void init_u( lapack_int size, lapack_complex_double *u );
static void init_c( lapack_int size, lapack_complex_double *c );
static void init_work( lapack_int size, double *work );
static int compare_zbdsqr( double *d, double *d_i, double *e, double *e_i,
                           lapack_complex_double *vt,
                           lapack_complex_double *vt_i,
                           lapack_complex_double *u, lapack_complex_double *u_i,
                           lapack_complex_double *c, lapack_complex_double *c_i,
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
    double *d = NULL, *d_i = NULL;
    double *e = NULL, *e_i = NULL;
    lapack_complex_double *vt = NULL, *vt_i = NULL;
    lapack_complex_double *u = NULL, *u_i = NULL;
    lapack_complex_double *c = NULL, *c_i = NULL;
    double *work = NULL, *work_i = NULL;
    double *d_save = NULL;
    double *e_save = NULL;
    lapack_complex_double *vt_save = NULL;
    lapack_complex_double *u_save = NULL;
    lapack_complex_double *c_save = NULL;
    lapack_complex_double *vt_r = NULL;
    lapack_complex_double *u_r = NULL;
    lapack_complex_double *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_zbdsqr( &uplo, &n, &ncvt, &nru, &ncc, &ldvt, &ldu, &ldc );
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
    d = (double *)LAPACKE_malloc( n * sizeof(double) );
    e = (double *)LAPACKE_malloc( n * sizeof(double) );
    vt = (lapack_complex_double *)
        LAPACKE_malloc( ldvt*ncvt * sizeof(lapack_complex_double) );
    u = (lapack_complex_double *)
        LAPACKE_malloc( ldu*n * sizeof(lapack_complex_double) );
    c = (lapack_complex_double *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_double) );
    work = (double *)LAPACKE_malloc( 4*n * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    d_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    e_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    vt_i = (lapack_complex_double *)
        LAPACKE_malloc( ldvt*ncvt * sizeof(lapack_complex_double) );
    u_i = (lapack_complex_double *)
        LAPACKE_malloc( ldu*n * sizeof(lapack_complex_double) );
    c_i = (lapack_complex_double *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_double) );
    work_i = (double *)LAPACKE_malloc( 4*n * sizeof(double) );

    /* Allocate memory for the backup arrays */
    d_save = (double *)LAPACKE_malloc( n * sizeof(double) );
    e_save = (double *)LAPACKE_malloc( n * sizeof(double) );
    vt_save = (lapack_complex_double *)
        LAPACKE_malloc( ldvt*ncvt * sizeof(lapack_complex_double) );
    u_save = (lapack_complex_double *)
        LAPACKE_malloc( ldu*n * sizeof(lapack_complex_double) );
    c_save = (lapack_complex_double *)
        LAPACKE_malloc( ldc*ncc * sizeof(lapack_complex_double) );

    /* Allocate memory for the row-major arrays */
    vt_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(ncvt+2) * sizeof(lapack_complex_double) );
    u_r = (lapack_complex_double *)
        LAPACKE_malloc( nru*(n+2) * sizeof(lapack_complex_double) );
    c_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(ncc+2) * sizeof(lapack_complex_double) );

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
    zbdsqr_( &uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc,
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
    info_i = LAPACKE_zbdsqr_work( LAPACK_COL_MAJOR, uplo_i, n_i, ncvt_i, nru_i,
                                  ncc_i, d_i, e_i, vt_i, ldvt_i, u_i, ldu_i,
                                  c_i, ldc_i, work_i );

    failed = compare_zbdsqr( d, d_i, e, e_i, vt, vt_i, u, u_i, c, c_i, info,
                             info_i, ldc, ldu, ldvt, n, ncc, ncvt, nru );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zbdsqr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zbdsqr\n" );
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
    info_i = LAPACKE_zbdsqr( LAPACK_COL_MAJOR, uplo_i, n_i, ncvt_i, nru_i,
                             ncc_i, d_i, e_i, vt_i, ldvt_i, u_i, ldu_i, c_i,
                             ldc_i );

    failed = compare_zbdsqr( d, d_i, e, e_i, vt, vt_i, u, u_i, c, c_i, info,
                             info_i, ldc, ldu, ldvt, n, ncc, ncvt, nru );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zbdsqr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zbdsqr\n" );
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
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, ncvt, vt_i, ldvt, vt_r,
                           ncvt+2 );
    }
    if( nru != 0 ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, nru, n, u_i, ldu, u_r, n+2 );
    }
    if( ncc != 0 ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, ncc, c_i, ldc, c_r, ncc+2 );
    }
    info_i = LAPACKE_zbdsqr_work( LAPACK_ROW_MAJOR, uplo_i, n_i, ncvt_i, nru_i,
                                  ncc_i, d_i, e_i, vt_r, ldvt_r, u_r, ldu_r,
                                  c_r, ldc_r, work_i );

    if( ncvt != 0 ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, ncvt, vt_r, ncvt+2, vt_i,
                           ldvt );
    }
    if( nru != 0 ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, nru, n, u_r, n+2, u_i, ldu );
    }
    if( ncc != 0 ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, ncc, c_r, ncc+2, c_i, ldc );
    }

    failed = compare_zbdsqr( d, d_i, e, e_i, vt, vt_i, u, u_i, c, c_i, info,
                             info_i, ldc, ldu, ldvt, n, ncc, ncvt, nru );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zbdsqr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zbdsqr\n" );
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
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, ncvt, vt_i, ldvt, vt_r,
                           ncvt+2 );
    }
    if( nru != 0 ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, nru, n, u_i, ldu, u_r, n+2 );
    }
    if( ncc != 0 ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, ncc, c_i, ldc, c_r, ncc+2 );
    }
    info_i = LAPACKE_zbdsqr( LAPACK_ROW_MAJOR, uplo_i, n_i, ncvt_i, nru_i,
                             ncc_i, d_i, e_i, vt_r, ldvt_r, u_r, ldu_r, c_r,
                             ldc_r );

    if( ncvt != 0 ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, ncvt, vt_r, ncvt+2, vt_i,
                           ldvt );
    }
    if( nru != 0 ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, nru, n, u_r, n+2, u_i, ldu );
    }
    if( ncc != 0 ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, ncc, c_r, ncc+2, c_i, ldc );
    }

    failed = compare_zbdsqr( d, d_i, e, e_i, vt, vt_i, u, u_i, c, c_i, info,
                             info_i, ldc, ldu, ldvt, n, ncc, ncvt, nru );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zbdsqr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zbdsqr\n" );
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

/* Auxiliary function: zbdsqr scalar parameters initialization */
static void init_scalars_zbdsqr( char *uplo, lapack_int *n, lapack_int *ncvt,
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

/* Auxiliary functions: zbdsqr array parameters initialization */
static void init_d( lapack_int size, double *d ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        d[i] = 0;
    }
    d[0] = -2.22551117723546850e+000;
    d[1] = 1.60455842542666320e+000;
    d[2] = -1.67305626481148150e+000;
}
static void init_e( lapack_int size, double *e ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        e[i] = 0;
    }
    e[0] = -1.63488453997962370e+000;
    e[1] = 9.91931214905034660e-001;
    e[2] = 0.00000000000000000e+000;
}
static void init_vt( lapack_int size, lapack_complex_double *vt ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vt[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    vt[0] = lapack_make_complex_double( -1.25813791844360030e-001,
                                        1.61760589514177240e-001 );
    vt[8] = lapack_make_complex_double( -2.24667485436357280e-001,
                                        3.86428074950534490e-001 );
    vt[16] = lapack_make_complex_double( 3.45987927571990260e-001,
                                         2.15680786018902980e-001 );
    vt[24] = lapack_make_complex_double( -7.09949253978889080e-001,
                                         -2.96561080775991640e-001 );
    vt[1] = lapack_make_complex_double( 4.14785092219841440e-001,
                                        1.79504972515068440e-001 );
    vt[9] = lapack_make_complex_double( 1.36782286410614380e-001,
                                        -3.97619520383010570e-001 );
    vt[17] = lapack_make_complex_double( 6.88479208824614890e-001,
                                         3.38638736715237810e-001 );
    vt[25] = lapack_make_complex_double( 1.66728569105729860e-001,
                                         -4.94204866963832150e-002 );
    vt[2] = lapack_make_complex_double( 4.57518505565114970e-001,
                                        -4.80657061196396170e-001 );
    vt[10] = lapack_make_complex_double( -2.73326699295941440e-001,
                                         4.98067456784206160e-001 );
    vt[18] = lapack_make_complex_double( -2.29670751731120840e-002,
                                         3.86052844703753340e-001 );
    vt[26] = lapack_make_complex_double( 1.73047045142453560e-001,
                                         2.39493977489214750e-001 );
}
static void init_u( lapack_int size, lapack_complex_double *u ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        u[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    u[0] = lapack_make_complex_double( 1.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    u[8] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    u[16] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    u[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    u[9] = lapack_make_complex_double( -5.02026435338565060e-001,
                                       -7.57519577948184250e-001 );
    u[17] = lapack_make_complex_double( -3.59686703621296070e-001,
                                        2.11563282458238110e-001 );
    u[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    u[10] = lapack_make_complex_double( -6.32136648022973980e-004,
                                        4.17292640301403450e-001 );
    u[18] = lapack_make_complex_double( -9.07393061823706050e-001,
                                        -5.00428227395494780e-002 );
}
static void init_c( lapack_int size, lapack_complex_double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to zbdsqr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zbdsqr( double *d, double *d_i, double *e, double *e_i,
                           lapack_complex_double *vt,
                           lapack_complex_double *vt_i,
                           lapack_complex_double *u, lapack_complex_double *u_i,
                           lapack_complex_double *c, lapack_complex_double *c_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int ldu, lapack_int ldvt, lapack_int n,
                           lapack_int ncc, lapack_int ncvt, lapack_int nru )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < n; i++ ) {
        failed += compare_doubles(d[i],d_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_doubles(e[i],e_i[i]);
    }
    if( ncvt != 0 ) {
        for( i = 0; i < ldvt*ncvt; i++ ) {
            failed += compare_complex_doubles(vt[i],vt_i[i]);
        }
    }
    if( nru != 0 ) {
        for( i = 0; i < ldu*n; i++ ) {
            failed += compare_complex_doubles(u[i],u_i[i]);
        }
    }
    if( ncc != 0 ) {
        for( i = 0; i < ldc*ncc; i++ ) {
            failed += compare_complex_doubles(c[i],c_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
