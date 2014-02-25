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
* ctrsyl_1 is the test program for the C interface to LAPACK
* routine ctrsyl
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

static void init_scalars_ctrsyl( char *trana, char *tranb, lapack_int *isgn,
                                 lapack_int *m, lapack_int *n, lapack_int *lda,
                                 lapack_int *ldb, lapack_int *ldc );
static void init_a( lapack_int size, lapack_complex_float *a );
static void init_b( lapack_int size, lapack_complex_float *b );
static void init_c( lapack_int size, lapack_complex_float *c );
static int compare_ctrsyl( lapack_complex_float *c, lapack_complex_float *c_i,
                           float scale, float scale_i, lapack_int info,
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
    float scale, scale_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *a = NULL, *a_i = NULL;
    lapack_complex_float *b = NULL, *b_i = NULL;
    lapack_complex_float *c = NULL, *c_i = NULL;
    lapack_complex_float *c_save = NULL;
    lapack_complex_float *a_r = NULL;
    lapack_complex_float *b_r = NULL;
    lapack_complex_float *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ctrsyl( &trana, &tranb, &isgn, &m, &n, &lda, &ldb, &ldc );
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
    a = (lapack_complex_float *)
        LAPACKE_malloc( lda*m * sizeof(lapack_complex_float) );
    b = (lapack_complex_float *)
        LAPACKE_malloc( ldb*n * sizeof(lapack_complex_float) );
    c = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_float *)
        LAPACKE_malloc( lda*m * sizeof(lapack_complex_float) );
    b_i = (lapack_complex_float *)
        LAPACKE_malloc( ldb*n * sizeof(lapack_complex_float) );
    c_i = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    c_save = (lapack_complex_float *)
        LAPACKE_malloc( ldc*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_float *)
        LAPACKE_malloc( m*(m+2) * sizeof(lapack_complex_float) );
    b_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );
    c_r = (lapack_complex_float *)
        LAPACKE_malloc( m*(n+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_a( lda*m, a );
    init_b( ldb*n, b );
    init_c( ldc*n, c );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldc*n; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    ctrsyl_( &trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, &scale,
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
    info_i = LAPACKE_ctrsyl_work( LAPACK_COL_MAJOR, trana_i, tranb_i, isgn_i,
                                  m_i, n_i, a_i, lda_i, b_i, ldb_i, c_i, ldc_i,
                                  &scale_i );

    failed = compare_ctrsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ctrsyl\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ctrsyl\n" );
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
    info_i = LAPACKE_ctrsyl( LAPACK_COL_MAJOR, trana_i, tranb_i, isgn_i, m_i,
                             n_i, a_i, lda_i, b_i, ldb_i, c_i, ldc_i,
                             &scale_i );

    failed = compare_ctrsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ctrsyl\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ctrsyl\n" );
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

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, m, a_i, lda, a_r, m+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, b_i, ldb, b_r, n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_ctrsyl_work( LAPACK_ROW_MAJOR, trana_i, tranb_i, isgn_i,
                                  m_i, n_i, a_r, lda_r, b_r, ldb_r, c_r, ldc_r,
                                  &scale_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_ctrsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ctrsyl\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ctrsyl\n" );
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
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, m, a_i, lda, a_r, m+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, b_i, ldb, b_r, n+2 );
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_ctrsyl( LAPACK_ROW_MAJOR, trana_i, tranb_i, isgn_i, m_i,
                             n_i, a_r, lda_r, b_r, ldb_r, c_r, ldc_r,
                             &scale_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_ctrsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ctrsyl\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ctrsyl\n" );
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

/* Auxiliary function: ctrsyl scalar parameters initialization */
static void init_scalars_ctrsyl( char *trana, char *tranb, lapack_int *isgn,
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

/* Auxiliary functions: ctrsyl array parameters initialization */
static void init_a( lapack_int size, lapack_complex_float *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    a[0] = lapack_make_complex_float( -6.000000000e+000, -7.000000000e+000 );
    a[8] = lapack_make_complex_float( 3.600000143e-001, -3.600000143e-001 );
    a[16] = lapack_make_complex_float( -1.899999976e-001, 4.799999893e-001 );
    a[24] = lapack_make_complex_float( 8.799999952e-001, -2.500000000e-001 );
    a[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[9] = lapack_make_complex_float( -5.000000000e+000, 2.000000000e+000 );
    a[17] = lapack_make_complex_float( -2.999999933e-002, -7.200000286e-001 );
    a[25] = lapack_make_complex_float( -2.300000042e-001, 1.299999952e-001 );
    a[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[10] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[18] = lapack_make_complex_float( 8.000000000e+000, -1.000000000e+000 );
    a[26] = lapack_make_complex_float( 9.399999976e-001, 5.299999714e-001 );
    a[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[19] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[27] = lapack_make_complex_float( 3.000000000e+000, -4.000000000e+000 );
}
static void init_b( lapack_int size, lapack_complex_float *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    b[0] = lapack_make_complex_float( 5.000000000e-001, -2.000000030e-001 );
    b[8] = lapack_make_complex_float( -2.899999917e-001, -1.599999964e-001 );
    b[16] = lapack_make_complex_float( -3.700000048e-001, 8.399999738e-001 );
    b[24] = lapack_make_complex_float( -5.500000119e-001, 7.300000191e-001 );
    b[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    b[9] = lapack_make_complex_float( -4.000000060e-001, 8.999999762e-001 );
    b[17] = lapack_make_complex_float( 5.999999866e-002, 2.199999988e-001 );
    b[25] = lapack_make_complex_float( -4.300000072e-001, 1.700000018e-001 );
    b[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    b[10] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    b[18] = lapack_make_complex_float( -8.999999762e-001, -1.000000015e-001 );
    b[26] = lapack_make_complex_float( -8.899999857e-001, -4.199999869e-001 );
    b[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    b[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    b[19] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    b[27] = lapack_make_complex_float( 3.000000119e-001, -6.999999881e-001 );
}
static void init_c( lapack_int size, lapack_complex_float *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    c[0] = lapack_make_complex_float( 6.299999952e-001, 3.499999940e-001 );
    c[8] = lapack_make_complex_float( 4.499999881e-001, -5.600000024e-001 );
    c[16] = lapack_make_complex_float( 7.999999821e-002, -1.400000006e-001 );
    c[24] = lapack_make_complex_float( -1.700000018e-001, -2.300000042e-001 );
    c[1] = lapack_make_complex_float( -1.700000018e-001, 9.000000358e-002 );
    c[9] = lapack_make_complex_float( -7.000000030e-002, -3.100000024e-001 );
    c[17] = lapack_make_complex_float( 2.700000107e-001, -5.400000215e-001 );
    c[25] = lapack_make_complex_float( 3.499999940e-001, 1.210000038e+000 );
    c[2] = lapack_make_complex_float( -9.300000072e-001, -4.399999976e-001 );
    c[10] = lapack_make_complex_float( -3.300000131e-001, -3.499999940e-001 );
    c[18] = lapack_make_complex_float( 4.099999964e-001, -2.999999933e-002 );
    c[26] = lapack_make_complex_float( 5.699999928e-001, 8.399999738e-001 );
    c[3] = lapack_make_complex_float( 5.400000215e-001, 2.500000000e-001 );
    c[11] = lapack_make_complex_float( -6.200000048e-001, -5.000000075e-002 );
    c[19] = lapack_make_complex_float( -5.199999809e-001, -1.299999952e-001 );
    c[27] = lapack_make_complex_float( 1.099999994e-001, -7.999999821e-002 );
}

/* Auxiliary function: C interface to ctrsyl results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ctrsyl( lapack_complex_float *c, lapack_complex_float *c_i,
                           float scale, float scale_i, lapack_int info,
                           lapack_int info_i, lapack_int ldc, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldc*n; i++ ) {
        failed += compare_complex_floats(c[i],c_i[i]);
    }
    failed += compare_floats(scale,scale_i);
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
