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
* ztrsyl_1 is the test program for the C interface to LAPACK
* routine ztrsyl
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

static void init_scalars_ztrsyl( char *trana, char *tranb, lapack_int *isgn,
                                 lapack_int *m, lapack_int *n, lapack_int *lda,
                                 lapack_int *ldb, lapack_int *ldc );
static void init_a( lapack_int size, lapack_complex_double *a );
static void init_b( lapack_int size, lapack_complex_double *b );
static void init_c( lapack_int size, lapack_complex_double *c );
static int compare_ztrsyl( lapack_complex_double *c, lapack_complex_double *c_i,
                           double scale, double scale_i, lapack_int info,
                           lapack_int info_i, lapack_int ldc, lapack_int n );

int main(void)
{
    /* Local scalars */
    char trana, trana_i;
    char tranb, tranb_i;
    lapack_int isgn, isgn_i;
    lapack_int m, m_i;
    lapack_int n, n_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    lapack_int ldb, ldb_i;
    lapack_int ldb_r;
    lapack_int ldc, ldc_i;
    lapack_int ldc_r;
    double scale, scale_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_double *a = NULL, *a_i = NULL;
    lapack_complex_double *b = NULL, *b_i = NULL;
    lapack_complex_double *c = NULL, *c_i = NULL;
    lapack_complex_double *c_save = NULL;
    lapack_complex_double *a_r = NULL;
    lapack_complex_double *b_r = NULL;
    lapack_complex_double *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ztrsyl( &trana, &tranb, &isgn, &m, &n, &lda, &ldb, &ldc );
    lda_r = m+2;
    ldb_r = n+2;
    ldc_r = n+2;
    trana_i = trana;
    tranb_i = tranb;
    isgn_i = isgn;
    m_i = m;
    n_i = n;
    lda_i = lda;
    ldb_i = ldb;
    ldc_i = ldc;

    /* Allocate memory for the LAPACK routine arrays */
    a = (lapack_complex_double *)
        LAPACKE_malloc( lda*m * sizeof(lapack_complex_double) );
    b = (lapack_complex_double *)
        LAPACKE_malloc( ldb*n * sizeof(lapack_complex_double) );
    c = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_double *)
        LAPACKE_malloc( lda*m * sizeof(lapack_complex_double) );
    b_i = (lapack_complex_double *)
        LAPACKE_malloc( ldb*n * sizeof(lapack_complex_double) );
    c_i = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the backup arrays */
    c_save = (lapack_complex_double *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_double *)
        LAPACKE_malloc( m*(m+2) * sizeof(lapack_complex_double) );
    b_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );
    c_r = (lapack_complex_double *)
        LAPACKE_malloc( m*(n+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_a( lda*m, a );
    init_b( ldb*n, b );
    init_c( ldc*n, c );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldc*n; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    ztrsyl_( &trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, &scale,
             &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*m; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < ldb*n; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    info_i = LAPACKE_ztrsyl_work( LAPACK_COL_MAJOR, trana_i, tranb_i, isgn_i,
                                  m_i, n_i, a_i, lda_i, b_i, ldb_i, c_i, ldc_i,
                                  &scale_i );

    failed = compare_ztrsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ztrsyl\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ztrsyl\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*m; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < ldb*n; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    info_i = LAPACKE_ztrsyl( LAPACK_COL_MAJOR, trana_i, tranb_i, isgn_i, m_i,
                             n_i, a_i, lda_i, b_i, ldb_i, c_i, ldc_i,
                             &scale_i );

    failed = compare_ztrsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ztrsyl\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ztrsyl\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*m; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < ldb*n; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, m, a_i, lda, a_r, m+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, b_i, ldb, b_r, n+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_ztrsyl_work( LAPACK_ROW_MAJOR, trana_i, tranb_i, isgn_i,
                                  m_i, n_i, a_r, lda_r, b_r, ldb_r, c_r, ldc_r,
                                  &scale_i );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_ztrsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ztrsyl\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ztrsyl\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*m; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < ldb*n; i++ ) {
        b_i[i] = b[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, m, a_i, lda, a_r, m+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, b_i, ldb, b_r, n+2 );
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_ztrsyl( LAPACK_ROW_MAJOR, trana_i, tranb_i, isgn_i, m_i,
                             n_i, a_r, lda_r, b_r, ldb_r, c_r, ldc_r,
                             &scale_i );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_ztrsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ztrsyl\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ztrsyl\n" );
    }

    /* Release memory */
    if( a != NULL ) {
        LAPACKE_free( a );
    }
    if( a_i != NULL ) {
        LAPACKE_free( a_i );
    }
    if( a_r != NULL ) {
        LAPACKE_free( a_r );
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

    return 0;
}

/* Auxiliary function: ztrsyl scalar parameters initialization */
static void init_scalars_ztrsyl( char *trana, char *tranb, lapack_int *isgn,
                                 lapack_int *m, lapack_int *n, lapack_int *lda,
                                 lapack_int *ldb, lapack_int *ldc )
{
    *trana = 'N';
    *tranb = 'N';
    *isgn = 1;
    *m = 4;
    *n = 4;
    *lda = 8;
    *ldb = 8;
    *ldc = 8;

    return;
}

/* Auxiliary functions: ztrsyl array parameters initialization */
static void init_a( lapack_int size, lapack_complex_double *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    a[0] = lapack_make_complex_double( -6.00000000000000000e+000,
                                       -7.00000000000000000e+000 );
    a[8] = lapack_make_complex_double( 3.59999999999999990e-001,
                                       -3.59999999999999990e-001 );
    a[16] = lapack_make_complex_double( -1.90000000000000000e-001,
                                        4.79999999999999980e-001 );
    a[24] = lapack_make_complex_double( 8.80000000000000000e-001,
                                        -2.50000000000000000e-001 );
    a[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    a[9] = lapack_make_complex_double( -5.00000000000000000e+000,
                                       2.00000000000000000e+000 );
    a[17] = lapack_make_complex_double( -2.99999999999999990e-002,
                                        -7.19999999999999970e-001 );
    a[25] = lapack_make_complex_double( -2.30000000000000010e-001,
                                        1.30000000000000000e-001 );
    a[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    a[10] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    a[18] = lapack_make_complex_double( 8.00000000000000000e+000,
                                        -1.00000000000000000e+000 );
    a[26] = lapack_make_complex_double( 9.39999999999999950e-001,
                                        5.30000000000000030e-001 );
    a[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    a[11] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    a[19] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    a[27] = lapack_make_complex_double( 3.00000000000000000e+000,
                                        -4.00000000000000000e+000 );
}
static void init_b( lapack_int size, lapack_complex_double *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    b[0] = lapack_make_complex_double( 5.00000000000000000e-001,
                                       -2.00000000000000010e-001 );
    b[8] = lapack_make_complex_double( -2.89999999999999980e-001,
                                       -1.60000000000000000e-001 );
    b[16] = lapack_make_complex_double( -3.70000000000000000e-001,
                                        8.39999999999999970e-001 );
    b[24] = lapack_make_complex_double( -5.50000000000000040e-001,
                                        7.29999999999999980e-001 );
    b[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    b[9] = lapack_make_complex_double( -4.00000000000000020e-001,
                                       9.00000000000000020e-001 );
    b[17] = lapack_make_complex_double( 5.99999999999999980e-002,
                                        2.20000000000000000e-001 );
    b[25] = lapack_make_complex_double( -4.29999999999999990e-001,
                                        1.70000000000000010e-001 );
    b[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    b[10] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    b[18] = lapack_make_complex_double( -9.00000000000000020e-001,
                                        -1.00000000000000010e-001 );
    b[26] = lapack_make_complex_double( -8.90000000000000010e-001,
                                        -4.19999999999999980e-001 );
    b[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    b[11] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    b[19] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    b[27] = lapack_make_complex_double( 2.99999999999999990e-001,
                                        -6.99999999999999960e-001 );
}
static void init_c( lapack_int size, lapack_complex_double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    c[0] = lapack_make_complex_double( 6.30000000000000000e-001,
                                       3.49999999999999980e-001 );
    c[8] = lapack_make_complex_double( 4.50000000000000010e-001,
                                       -5.60000000000000050e-001 );
    c[16] = lapack_make_complex_double( 8.00000000000000020e-002,
                                        -1.40000000000000010e-001 );
    c[24] = lapack_make_complex_double( -1.70000000000000010e-001,
                                        -2.30000000000000010e-001 );
    c[1] = lapack_make_complex_double( -1.70000000000000010e-001,
                                       8.99999999999999970e-002 );
    c[9] = lapack_make_complex_double( -7.00000000000000070e-002,
                                       -3.10000000000000000e-001 );
    c[17] = lapack_make_complex_double( 2.70000000000000020e-001,
                                        -5.40000000000000040e-001 );
    c[25] = lapack_make_complex_double( 3.49999999999999980e-001,
                                        1.21000000000000000e+000 );
    c[2] = lapack_make_complex_double( -9.30000000000000050e-001,
                                       -4.40000000000000000e-001 );
    c[10] = lapack_make_complex_double( -3.30000000000000020e-001,
                                        -3.49999999999999980e-001 );
    c[18] = lapack_make_complex_double( 4.09999999999999980e-001,
                                        -2.99999999999999990e-002 );
    c[26] = lapack_make_complex_double( 5.69999999999999950e-001,
                                        8.39999999999999970e-001 );
    c[3] = lapack_make_complex_double( 5.40000000000000040e-001,
                                       2.50000000000000000e-001 );
    c[11] = lapack_make_complex_double( -6.20000000000000000e-001,
                                        -5.00000000000000030e-002 );
    c[19] = lapack_make_complex_double( -5.20000000000000020e-001,
                                        -1.30000000000000000e-001 );
    c[27] = lapack_make_complex_double( 1.10000000000000000e-001,
                                        -8.00000000000000020e-002 );
}

/* Auxiliary function: C interface to ztrsyl results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ztrsyl( lapack_complex_double *c, lapack_complex_double *c_i,
                           double scale, double scale_i, lapack_int info,
                           lapack_int info_i, lapack_int ldc, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldc*n; i++ ) {
        failed += compare_complex_doubles(c[i],c_i[i]);
    }
    failed += compare_doubles(scale,scale_i);
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
