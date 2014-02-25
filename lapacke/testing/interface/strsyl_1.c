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
* strsyl_1 is the test program for the C interface to LAPACK
* routine strsyl
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

static void init_scalars_strsyl( char *trana, char *tranb, lapack_int *isgn,
                                 lapack_int *m, lapack_int *n, lapack_int *lda,
                                 lapack_int *ldb, lapack_int *ldc );
static void init_a( lapack_int size, float *a );
static void init_b( lapack_int size, float *b );
static void init_c( lapack_int size, float *c );
static int compare_strsyl( float *c, float *c_i, float scale, float scale_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int n );

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
    float *a = NULL, *a_i = NULL;
    float *b = NULL, *b_i = NULL;
    float *c = NULL, *c_i = NULL;
    float *c_save = NULL;
    float *a_r = NULL;
    float *b_r = NULL;
    float *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_strsyl( &trana, &tranb, &isgn, &m, &n, &lda, &ldb, &ldc );
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
    a = (float *)LAPACKE_malloc( lda*m * sizeof(float) );
    b = (float *)LAPACKE_malloc( ldb*n * sizeof(float) );
    c = (float *)LAPACKE_malloc( ldc*n * sizeof(float) );

    /* Allocate memory for the C interface function arrays */
    a_i = (float *)LAPACKE_malloc( lda*m * sizeof(float) );
    b_i = (float *)LAPACKE_malloc( ldb*n * sizeof(float) );
    c_i = (float *)LAPACKE_malloc( ldc*n * sizeof(float) );

    /* Allocate memory for the backup arrays */
    c_save = (float *)LAPACKE_malloc( ldc*n * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    a_r = (float *)LAPACKE_malloc( m*(m+2) * sizeof(float) );
    b_r = (float *)LAPACKE_malloc( n*(n+2) * sizeof(float) );
    c_r = (float *)LAPACKE_malloc( m*(n+2) * sizeof(float) );

    /* Initialize input arrays */
    init_a( lda*m, a );
    init_b( ldb*n, b );
    init_c( ldc*n, c );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldc*n; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    strsyl_( &trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, &scale,
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
    info_i = LAPACKE_strsyl_work( LAPACK_COL_MAJOR, trana_i, tranb_i, isgn_i,
                                  m_i, n_i, a_i, lda_i, b_i, ldb_i, c_i, ldc_i,
                                  &scale_i );

    failed = compare_strsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to strsyl\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to strsyl\n" );
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
    info_i = LAPACKE_strsyl( LAPACK_COL_MAJOR, trana_i, tranb_i, isgn_i, m_i,
                             n_i, a_i, lda_i, b_i, ldb_i, c_i, ldc_i,
                             &scale_i );

    failed = compare_strsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to strsyl\n" );
    } else {
        printf( "FAILED: column-major high-level interface to strsyl\n" );
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

    LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, m, a_i, lda, a_r, m+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, b_i, ldb, b_r, n+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_strsyl_work( LAPACK_ROW_MAJOR, trana_i, tranb_i, isgn_i,
                                  m_i, n_i, a_r, lda_r, b_r, ldb_r, c_r, ldc_r,
                                  &scale_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_strsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to strsyl\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to strsyl\n" );
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
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, m, a_i, lda, a_r, m+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, b_i, ldb, b_r, n+2 );
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_strsyl( LAPACK_ROW_MAJOR, trana_i, tranb_i, isgn_i, m_i,
                             n_i, a_r, lda_r, b_r, ldb_r, c_r, ldc_r,
                             &scale_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_strsyl( c, c_i, scale, scale_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to strsyl\n" );
    } else {
        printf( "FAILED: row-major high-level interface to strsyl\n" );
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

/* Auxiliary function: strsyl scalar parameters initialization */
static void init_scalars_strsyl( char *trana, char *tranb, lapack_int *isgn,
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

/* Auxiliary functions: strsyl array parameters initialization */
static void init_a( lapack_int size, float *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = 0;
    }
    a[0] = 1.000000015e-001;  /* a[0,0] */
    a[8] = 5.000000000e-001;  /* a[0,1] */
    a[16] = 6.800000072e-001;  /* a[0,2] */
    a[24] = -2.099999934e-001;  /* a[0,3] */
    a[1] = -5.000000000e-001;  /* a[1,0] */
    a[9] = 1.000000015e-001;  /* a[1,1] */
    a[17] = -2.399999946e-001;  /* a[1,2] */
    a[25] = 6.700000167e-001;  /* a[1,3] */
    a[2] = 0.000000000e+000;  /* a[2,0] */
    a[10] = 0.000000000e+000;  /* a[2,1] */
    a[18] = 1.899999976e-001;  /* a[2,2] */
    a[26] = -3.499999940e-001;  /* a[2,3] */
    a[3] = 0.000000000e+000;  /* a[3,0] */
    a[11] = 0.000000000e+000;  /* a[3,1] */
    a[19] = 0.000000000e+000;  /* a[3,2] */
    a[27] = -7.200000286e-001;  /* a[3,3] */
}
static void init_b( lapack_int size, float *b ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        b[i] = 0;
    }
    b[0] = -9.900000095e-001;  /* b[0,0] */
    b[8] = -1.700000018e-001;  /* b[0,1] */
    b[16] = 3.899999857e-001;  /* b[0,2] */
    b[24] = 5.799999833e-001;  /* b[0,3] */
    b[1] = 0.000000000e+000;  /* b[1,0] */
    b[9] = 4.799999893e-001;  /* b[1,1] */
    b[17] = -8.399999738e-001;  /* b[1,2] */
    b[25] = -1.500000060e-001;  /* b[1,3] */
    b[2] = 0.000000000e+000;  /* b[2,0] */
    b[10] = 0.000000000e+000;  /* b[2,1] */
    b[18] = 7.500000000e-001;  /* b[2,2] */
    b[26] = 2.500000000e-001;  /* b[2,3] */
    b[3] = 0.000000000e+000;  /* b[3,0] */
    b[11] = 0.000000000e+000;  /* b[3,1] */
    b[19] = -2.500000000e-001;  /* b[3,2] */
    b[27] = 7.500000000e-001;  /* b[3,3] */
}
static void init_c( lapack_int size, float *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = 0;
    }
    c[0] = 6.299999952e-001;  /* c[0,0] */
    c[8] = -5.600000024e-001;  /* c[0,1] */
    c[16] = 7.999999821e-002;  /* c[0,2] */
    c[24] = -2.300000042e-001;  /* c[0,3] */
    c[1] = -4.499999881e-001;  /* c[1,0] */
    c[9] = -3.100000024e-001;  /* c[1,1] */
    c[17] = 2.700000107e-001;  /* c[1,2] */
    c[25] = 1.210000038e+000;  /* c[1,3] */
    c[2] = 2.000000030e-001;  /* c[2,0] */
    c[10] = -3.499999940e-001;  /* c[2,1] */
    c[18] = 4.099999964e-001;  /* c[2,2] */
    c[26] = 8.399999738e-001;  /* c[2,3] */
    c[3] = 4.900000095e-001;  /* c[3,0] */
    c[11] = -5.000000075e-002;  /* c[3,1] */
    c[19] = -5.199999809e-001;  /* c[3,2] */
    c[27] = -7.999999821e-002;  /* c[3,3] */
}

/* Auxiliary function: C interface to strsyl results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_strsyl( float *c, float *c_i, float scale, float scale_i,
                           lapack_int info, lapack_int info_i, lapack_int ldc,
                           lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldc*n; i++ ) {
        failed += compare_floats(c[i],c_i[i]);
    }
    failed += compare_floats(scale,scale_i);
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
