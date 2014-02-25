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
* ctrtri_1 is the test program for the C interface to LAPACK
* routine ctrtri
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

static void init_scalars_ctrtri( char *uplo, char *diag, lapack_int *n,
                                 lapack_int *lda );
static void init_a( lapack_int size, lapack_complex_float *a );
static int compare_ctrtri( lapack_complex_float *a, lapack_complex_float *a_i,
                           lapack_int info, lapack_int info_i, lapack_int lda,
                           lapack_int n );

int main(void)
{
    /* Local scalars */
    char uplo, uplo_i;
    char diag, diag_i;
    lapack_int n, n_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *a = NULL, *a_i = NULL;
    lapack_complex_float *a_save = NULL;
    lapack_complex_float *a_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ctrtri( &uplo, &diag, &n, &lda );
    lda_r = n+2;
    uplo_i = uplo;
    diag_i = diag;
    n_i = n;
    lda_i = lda;

    /* Allocate memory for the LAPACK routine arrays */
    a = (lapack_complex_float *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_float *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    a_save = (lapack_complex_float *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_a( lda*n, a );

    /* Backup the ouptut arrays */
    for( i = 0; i < lda*n; i++ ) {
        a_save[i] = a[i];
    }

    /* Call the LAPACK routine */
    ctrtri_( &uplo, &diag, &n, a, &lda, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    info_i = LAPACKE_ctrtri_work( LAPACK_COL_MAJOR, uplo_i, diag_i, n_i, a_i,
                                  lda_i );

    failed = compare_ctrtri( a, a_i, info, info_i, lda, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ctrtri\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ctrtri\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    info_i = LAPACKE_ctrtri( LAPACK_COL_MAJOR, uplo_i, diag_i, n_i, a_i,
                             lda_i );

    failed = compare_ctrtri( a, a_i, info, info_i, lda, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ctrtri\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ctrtri\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_ctrtri_work( LAPACK_ROW_MAJOR, uplo_i, diag_i, n_i, a_r,
                                  lda_r );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, a_r, n+2, a_i, lda );

    failed = compare_ctrtri( a, a_i, info, info_i, lda, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ctrtri\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ctrtri\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_ctrtri( LAPACK_ROW_MAJOR, uplo_i, diag_i, n_i, a_r,
                             lda_r );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, a_r, n+2, a_i, lda );

    failed = compare_ctrtri( a, a_i, info, info_i, lda, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ctrtri\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ctrtri\n" );
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
    if( a_save != NULL ) {
        LAPACKE_free( a_save );
    }

    return 0;
}

/* Auxiliary function: ctrtri scalar parameters initialization */
static void init_scalars_ctrtri( char *uplo, char *diag, lapack_int *n,
                                 lapack_int *lda )
{
    *uplo = 'L';
    *diag = 'N';
    *n = 4;
    *lda = 8;

    return;
}

/* Auxiliary functions: ctrtri array parameters initialization */
static void init_a( lapack_int size, lapack_complex_float *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    a[0] = lapack_make_complex_float( 4.780000210e+000, 4.559999943e+000 );
    a[8] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[16] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[24] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[1] = lapack_make_complex_float( 2.000000000e+000, -3.000000119e-001 );
    a[9] = lapack_make_complex_float( -4.110000134e+000, 1.250000000e+000 );
    a[17] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[25] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[2] = lapack_make_complex_float( 2.890000105e+000, -1.340000033e+000 );
    a[10] = lapack_make_complex_float( 2.359999895e+000, -4.250000000e+000 );
    a[18] = lapack_make_complex_float( 4.150000095e+000, 8.000000119e-001 );
    a[26] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    a[3] = lapack_make_complex_float( -1.889999986e+000, 1.149999976e+000 );
    a[11] = lapack_make_complex_float( 3.999999911e-002, -3.690000057e+000 );
    a[19] = lapack_make_complex_float( -1.999999955e-002, 4.600000083e-001 );
    a[27] = lapack_make_complex_float( 3.300000131e-001, -2.599999905e-001 );
}

/* Auxiliary function: C interface to ctrtri results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ctrtri( lapack_complex_float *a, lapack_complex_float *a_i,
                           lapack_int info, lapack_int info_i, lapack_int lda,
                           lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < lda*n; i++ ) {
        failed += compare_complex_floats(a[i],a_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
