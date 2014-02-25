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
* zbdsqr_1 is the test program for the C interface to LAPACK
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
    *uplo = 'U';
    *n = 4;
    *ncvt = 4;
    *nru = 6;
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
    d[0] = -3.08700502105195800e+000;
    d[1] = 2.06603927667906850e+000;
    d[2] = 1.87312889112571160e+000;
    d[3] = 2.00218286620699090e+000;
}
static void init_e( lapack_int size, double *e ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        e[i] = 0;
    }
    e[0] = 2.11257100745583950e+000;
    e[1] = 1.26281010665522400e+000;
    e[2] = -1.61263387280039260e+000;
    e[3] = 0.00000000000000000e+000;
}
static void init_vt( lapack_int size, lapack_complex_double *vt ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vt[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    vt[0] = lapack_make_complex_double( 1.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vt[8] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vt[16] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vt[24] = lapack_make_complex_double( 0.00000000000000000e+000,
                                         0.00000000000000000e+000 );
    vt[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vt[9] = lapack_make_complex_double( -2.31234531617601880e-001,
                                        -5.40426381454295000e-001 );
    vt[17] = lapack_make_complex_double( 1.78624075518133200e-001,
                                         -5.88728024331977130e-001 );
    vt[25] = lapack_make_complex_double( -4.04798435024739940e-001,
                                         -3.34814721344187110e-001 );
    vt[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vt[10] = lapack_make_complex_double( 4.96225984754380510e-001,
                                         2.64807608556462020e-001 );
    vt[18] = lapack_make_complex_double( -3.55623105209443960e-001,
                                         -6.88816765203852820e-001 );
    vt[26] = lapack_make_complex_double( 2.75655244497366060e-001,
                                         -8.19424169862172370e-002 );
    vt[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    vt[11] = lapack_make_complex_double( 2.91779089728414760e-001,
                                         5.02962804706486800e-001 );
    vt[19] = lapack_make_complex_double( 4.85070854407798900e-002,
                                         1.34920297541892970e-001 );
    vt[27] = lapack_make_complex_double( -7.15581919558924920e-001,
                                         -3.59554546978138210e-001 );
}
static void init_u( lapack_int size, lapack_complex_double *u ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        u[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    u[0] = lapack_make_complex_double( -3.10981029656006490e-001,
                                       2.62390243772255500e-001 );
    u[8] = lapack_make_complex_double( -6.52100156652580590e-001,
                                       -5.53232930919111300e-001 );
    u[16] = lapack_make_complex_double( -4.26652277391470560e-002,
                                        -3.60551875913458290e-002 );
    u[24] = lapack_make_complex_double( 2.63434291736194760e-001,
                                        7.41134830329821800e-002 );
    u[1] = lapack_make_complex_double( 3.17459801107173260e-001,
                                       -6.41398373665513440e-001 );
    u[9] = lapack_make_complex_double( -3.48794491065693380e-001,
                                       -7.20567339846693860e-002 );
    u[17] = lapack_make_complex_double( -2.28738568471465590e-001,
                                        -6.90884157766049650e-003 );
    u[25] = lapack_make_complex_double( -1.10140229907637420e-001,
                                        3.26180370964031520e-002 );
    u[2] = lapack_make_complex_double( -2.00841914986170850e-001,
                                       1.49011743376836450e-001 );
    u[10] = lapack_make_complex_double( 3.10251489958585980e-001,
                                        -2.30309142448547440e-002 );
    u[18] = lapack_make_complex_double( -1.85460167625633830e-001,
                                        1.81728049653940960e-001 );
    u[26] = lapack_make_complex_double( 2.95608864759607740e-001,
                                        -5.64781944942283730e-001 );
    u[3] = lapack_make_complex_double( 1.19857271846585860e-001,
                                       -1.23096657572169240e-001 );
    u[11] = lapack_make_complex_double( 4.59162169100148230e-003,
                                        4.64117081960858240e-004 );
    u[19] = lapack_make_complex_double( 3.30466205984792210e-001,
                                        -4.82071149828388450e-001 );
    u[27] = lapack_make_complex_double( 6.74932088473138290e-002,
                                        -3.46397021614317290e-001 );
    u[4] = lapack_make_complex_double( -2.68869015223422270e-001,
                                       -1.65208672004753480e-001 );
    u[12] = lapack_make_complex_double( -1.79380225742158510e-001,
                                        5.86074392732159330e-002 );
    u[20] = lapack_make_complex_double( 5.23537299523901270e-001,
                                        2.57977784560615260e-001 );
    u[28] = lapack_make_complex_double( -3.92695031623721160e-001,
                                        -1.44970740631150650e-001 );
    u[5] = lapack_make_complex_double( -3.49853658363007300e-001,
                                       9.07028003163352360e-002 );
    u[13] = lapack_make_complex_double( -8.28849784427782010e-002,
                                        5.05886798225579360e-002 );
    u[21] = lapack_make_complex_double( -3.20239252692931740e-001,
                                        -3.03796892849298670e-001 );
    u[29] = lapack_make_complex_double( -3.17357323819129010e-001,
                                        -3.24135324237070420e-001 );
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
