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
* chpgst_1 is the test program for the C interface to LAPACK
* routine chpgst
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

static void init_scalars_chpgst( lapack_int *itype, char *uplo, lapack_int *n );
static void init_ap( lapack_int size, lapack_complex_float *ap );
static void init_bp( lapack_int size, lapack_complex_float *bp );
static int compare_chpgst( lapack_complex_float *ap, lapack_complex_float *ap_i,
                           lapack_int info, lapack_int info_i, lapack_int n );

int main(void)
{
    /* Local scalars */
    lapack_int itype, itype_i;
    char uplo, uplo_i;
    lapack_int n, n_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *ap = NULL, *ap_i = NULL;
    lapack_complex_float *bp = NULL, *bp_i = NULL;
    lapack_complex_float *ap_save = NULL;
    lapack_complex_float *ap_r = NULL;
    lapack_complex_float *bp_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_chpgst( &itype, &uplo, &n );
    itype_i = itype;
    uplo_i = uplo;
    n_i = n;

    /* Allocate memory for the LAPACK routine arrays */
    ap = (lapack_complex_float *)
        LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(lapack_complex_float) );
    bp = (lapack_complex_float *)
        LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    ap_i = (lapack_complex_float *)
        LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(lapack_complex_float) );
    bp_i = (lapack_complex_float *)
        LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    ap_save = (lapack_complex_float *)
        LAPACKE_malloc( ((n*(n+1)/2)) * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    ap_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+1)/2 * sizeof(lapack_complex_float) );
    bp_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+1)/2 * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_ap( (n*(n+1)/2), ap );
    init_bp( (n*(n+1)/2), bp );

    /* Backup the ouptut arrays */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_save[i] = ap[i];
    }

    /* Call the LAPACK routine */
    chpgst_( &itype, &uplo, &n, ap, bp, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        bp_i[i] = bp[i];
    }
    info_i = LAPACKE_chpgst_work( LAPACK_COL_MAJOR, itype_i, uplo_i, n_i, ap_i,
                                  bp_i );

    failed = compare_chpgst( ap, ap_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to chpgst\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to chpgst\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        bp_i[i] = bp[i];
    }
    info_i = LAPACKE_chpgst( LAPACK_COL_MAJOR, itype_i, uplo_i, n_i, ap_i,
                             bp_i );

    failed = compare_chpgst( ap, ap_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to chpgst\n" );
    } else {
        printf( "FAILED: column-major high-level interface to chpgst\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        bp_i[i] = bp[i];
    }

    LAPACKE_cpp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    LAPACKE_cpp_trans( LAPACK_COL_MAJOR, uplo, n, bp_i, bp_r );
    info_i = LAPACKE_chpgst_work( LAPACK_ROW_MAJOR, itype_i, uplo_i, n_i, ap_r,
                                  bp_r );

    LAPACKE_cpp_trans( LAPACK_ROW_MAJOR, uplo, n, ap_r, ap_i );

    failed = compare_chpgst( ap, ap_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to chpgst\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to chpgst\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        ap_i[i] = ap_save[i];
    }
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        bp_i[i] = bp[i];
    }

    /* Init row_major arrays */
    LAPACKE_cpp_trans( LAPACK_COL_MAJOR, uplo, n, ap_i, ap_r );
    LAPACKE_cpp_trans( LAPACK_COL_MAJOR, uplo, n, bp_i, bp_r );
    info_i = LAPACKE_chpgst( LAPACK_ROW_MAJOR, itype_i, uplo_i, n_i, ap_r,
                             bp_r );

    LAPACKE_cpp_trans( LAPACK_ROW_MAJOR, uplo, n, ap_r, ap_i );

    failed = compare_chpgst( ap, ap_i, info, info_i, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to chpgst\n" );
    } else {
        printf( "FAILED: row-major high-level interface to chpgst\n" );
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
    if( bp != NULL ) {
        LAPACKE_free( bp );
    }
    if( bp_i != NULL ) {
        LAPACKE_free( bp_i );
    }
    if( bp_r != NULL ) {
        LAPACKE_free( bp_r );
    }

    return 0;
}

/* Auxiliary function: chpgst scalar parameters initialization */
static void init_scalars_chpgst( lapack_int *itype, char *uplo, lapack_int *n )
{
    *itype = 1;
    *uplo = 'L';
    *n = 4;

    return;
}

/* Auxiliary functions: chpgst array parameters initialization */
static void init_ap( lapack_int size, lapack_complex_float *ap ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ap[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    ap[0] = lapack_make_complex_float( -7.360000134e+000, 0.000000000e+000 );
    ap[1] = lapack_make_complex_float( 7.699999809e-001, 4.300000072e-001 );
    ap[2] = lapack_make_complex_float( -6.399999857e-001, 9.200000167e-001 );
    ap[3] = lapack_make_complex_float( 3.009999990e+000, 6.969999790e+000 );
    ap[4] = lapack_make_complex_float( 3.490000010e+000, 0.000000000e+000 );
    ap[5] = lapack_make_complex_float( 2.190000057e+000, -4.449999809e+000 );
    ap[6] = lapack_make_complex_float( 1.899999976e+000, -3.730000019e+000 );
    ap[7] = lapack_make_complex_float( 1.199999973e-001, 0.000000000e+000 );
    ap[8] = lapack_make_complex_float( 2.880000114e+000, 3.170000076e+000 );
    ap[9] = lapack_make_complex_float( -2.539999962e+000, 0.000000000e+000 );
}
static void init_bp( lapack_int size, lapack_complex_float *bp ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        bp[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    bp[0] = lapack_make_complex_float( 1.797220111e+000, 0.000000000e+000 );
    bp[1] = lapack_make_complex_float( 8.401864767e-001, 1.068316579e+000 );
    bp[2] = lapack_make_complex_float( 1.057188272e+000, -4.673885107e-001 );
    bp[3] = lapack_make_complex_float( 2.336942554e-001, -1.391037226e+000 );
    bp[4] = lapack_make_complex_float( 1.316353440e+000, 0.000000000e+000 );
    bp[5] = lapack_make_complex_float( -4.701749682e-001, 3.130658865e-001 );
    bp[6] = lapack_make_complex_float( 8.335255831e-002, 3.676066920e-002 );
    bp[7] = lapack_make_complex_float( 1.560392976e+000, 0.000000000e+000 );
    bp[8] = lapack_make_complex_float( 9.359616637e-001, 9.899691939e-001 );
    bp[9] = lapack_make_complex_float( 6.603332758e-001, 0.000000000e+000 );
}

/* Auxiliary function: C interface to chpgst results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_chpgst( lapack_complex_float *ap, lapack_complex_float *ap_i,
                           lapack_int info, lapack_int info_i, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < (n*(n+1)/2); i++ ) {
        failed += compare_complex_floats(ap[i],ap_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
