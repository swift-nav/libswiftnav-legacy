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
* cgecon_1 is the test program for the C interface to LAPACK
* routine cgecon
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

static void init_scalars_cgecon( char *norm, lapack_int *n, lapack_int *lda,
                                 float *anorm );
static void init_a( lapack_int size, lapack_complex_float *a );
static void init_work( lapack_int size, lapack_complex_float *work );
static void init_rwork( lapack_int size, float *rwork );
static int compare_cgecon( float rcond, float rcond_i, lapack_int info,
                           lapack_int info_i );

int main(void)
{
    /* Local scalars */
    char norm, norm_i;
    lapack_int n, n_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    float anorm, anorm_i;
    float rcond, rcond_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *a = NULL, *a_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    float *rwork = NULL, *rwork_i = NULL;
    lapack_complex_float *a_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_cgecon( &norm, &n, &lda, &anorm );
    lda_r = n+2;
    norm_i = norm;
    n_i = n;
    lda_i = lda;
    anorm_i = anorm;

    /* Allocate memory for the LAPACK routine arrays */
    a = (lapack_complex_float *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_float) );
    rwork = (float *)LAPACKE_malloc( 2*n * sizeof(float) );

    /* Allocate memory for the C interface function arrays */
    a_i = (lapack_complex_float *)
        LAPACKE_malloc( lda*n * sizeof(lapack_complex_float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( 2*n * sizeof(lapack_complex_float) );
    rwork_i = (float *)LAPACKE_malloc( 2*n * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    a_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_a( lda*n, a );
    init_work( 2*n, work );
    init_rwork( 2*n, rwork );

    /* Call the LAPACK routine */
    cgecon_( &norm, &n, a, &lda, &anorm, &rcond, work, rwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_cgecon_work( LAPACK_COL_MAJOR, norm_i, n_i, a_i, lda_i,
                                  anorm_i, &rcond_i, work_i, rwork_i );

    failed = compare_cgecon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to cgecon\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to cgecon\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        rwork_i[i] = rwork[i];
    }
    info_i = LAPACKE_cgecon( LAPACK_COL_MAJOR, norm_i, n_i, a_i, lda_i, anorm_i,
                             &rcond_i );

    failed = compare_cgecon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to cgecon\n" );
    } else {
        printf( "FAILED: column-major high-level interface to cgecon\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        rwork_i[i] = rwork[i];
    }

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_cgecon_work( LAPACK_ROW_MAJOR, norm_i, n_i, a_r, lda_r,
                                  anorm_i, &rcond_i, work_i, rwork_i );

    failed = compare_cgecon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to cgecon\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to cgecon\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < 2*n; i++ ) {
        rwork_i[i] = rwork[i];
    }

    /* Init row_major arrays */
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_cgecon( LAPACK_ROW_MAJOR, norm_i, n_i, a_r, lda_r, anorm_i,
                             &rcond_i );

    failed = compare_cgecon( rcond, rcond_i, info, info_i );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to cgecon\n" );
    } else {
        printf( "FAILED: row-major high-level interface to cgecon\n" );
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
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }
    if( rwork != NULL ) {
        LAPACKE_free( rwork );
    }
    if( rwork_i != NULL ) {
        LAPACKE_free( rwork_i );
    }

    return 0;
}

/* Auxiliary function: cgecon scalar parameters initialization */
static void init_scalars_cgecon( char *norm, lapack_int *n, lapack_int *lda,
                                 float *anorm )
{
    *norm = '1';
    *n = 4;
    *lda = 8;
    *anorm = 1.288382244e+001;

    return;
}

/* Auxiliary functions: cgecon array parameters initialization */
static void init_a( lapack_int size, lapack_complex_float *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    a[0] = lapack_make_complex_float( -3.289999962e+000, -2.390000105e+000 );
    a[8] = lapack_make_complex_float( -1.909999967e+000, 4.420000076e+000 );
    a[16] = lapack_make_complex_float( -1.400000006e-001, -1.350000024e+000 );
    a[24] = lapack_make_complex_float( 1.720000029e+000, 1.350000024e+000 );
    a[1] = lapack_make_complex_float( 2.376120239e-001, 2.559596300e-001 );
    a[9] = lapack_make_complex_float( 4.895180702e+000, -7.113622427e-001 );
    a[17] = lapack_make_complex_float( -4.622798264e-001, 1.696610570e+000 );
    a[25] = lapack_make_complex_float( 1.226852775e+000, 6.189731956e-001 );
    a[2] = lapack_make_complex_float( -1.019521058e-001, -7.010135055e-001 );
    a[10] = lapack_make_complex_float( -6.691496968e-001, 3.688699007e-001 );
    a[18] = lapack_make_complex_float( -5.141410828e+000, -1.129969597e+000 );
    a[26] = lapack_make_complex_float( 9.982581139e-001, 3.850151896e-001 );
    a[3] = lapack_make_complex_float( -5.358546376e-001, 2.707273066e-001 );
    a[11] = lapack_make_complex_float( -2.040217221e-001, 8.601181507e-001 );
    a[19] = lapack_make_complex_float( 8.233040571e-003, 1.210636795e-001 );
    a[27] = lapack_make_complex_float( 1.482392550e-001, -1.252242178e-001 );
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}
static void init_rwork( lapack_int size, float *rwork ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        rwork[i] = 0;
    }
}

/* Auxiliary function: C interface to cgecon results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_cgecon( float rcond, float rcond_i, lapack_int info,
                           lapack_int info_i )
{
    int failed = 0;
    failed += compare_floats(rcond,rcond_i);
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
