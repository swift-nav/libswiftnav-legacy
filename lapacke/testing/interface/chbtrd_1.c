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
* chbtrd_1 is the test program for the C interface to LAPACK
* routine chbtrd
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

static void init_scalars_chbtrd( char *vect, char *uplo, lapack_int *n,
                                 lapack_int *kd, lapack_int *ldab,
                                 lapack_int *ldq );
static void init_ab( lapack_int size, lapack_complex_float *ab );
static void init_d( lapack_int size, float *d );
static void init_e( lapack_int size, float *e );
static void init_q( lapack_int size, lapack_complex_float *q );
static void init_work( lapack_int size, lapack_complex_float *work );
static int compare_chbtrd( lapack_complex_float *ab, lapack_complex_float *ab_i,
                           float *d, float *d_i, float *e, float *e_i,
                           lapack_complex_float *q, lapack_complex_float *q_i,
                           lapack_int info, lapack_int info_i, lapack_int ldab,
                           lapack_int ldq, lapack_int n, char vect );

int main(void)
{
    /* Local scalars */
    char vect, vect_i;
    char uplo, uplo_i;
    lapack_int n, n_i;
    lapack_int kd, kd_i;
    lapack_int ldab, ldab_i;
    lapack_int ldab_r;
    lapack_int ldq, ldq_i;
    lapack_int ldq_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *ab = NULL, *ab_i = NULL;
    float *d = NULL, *d_i = NULL;
    float *e = NULL, *e_i = NULL;
    lapack_complex_float *q = NULL, *q_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    lapack_complex_float *ab_save = NULL;
    float *d_save = NULL;
    float *e_save = NULL;
    lapack_complex_float *q_save = NULL;
    lapack_complex_float *ab_r = NULL;
    lapack_complex_float *q_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_chbtrd( &vect, &uplo, &n, &kd, &ldab, &ldq );
    ldab_r = n+2;
    ldq_r = n+2;
    vect_i = vect;
    uplo_i = uplo;
    n_i = n;
    kd_i = kd;
    ldab_i = ldab;
    ldq_i = ldq;

    /* Allocate memory for the LAPACK routine arrays */
    ab = (lapack_complex_float *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_float) );
    d = (float *)LAPACKE_malloc( n * sizeof(float) );
    e = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );
    q = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    ab_i = (lapack_complex_float *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_float) );
    d_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    e_i = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );
    q_i = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    ab_save = (lapack_complex_float *)
        LAPACKE_malloc( ldab*n * sizeof(lapack_complex_float) );
    d_save = (float *)LAPACKE_malloc( n * sizeof(float) );
    e_save = (float *)LAPACKE_malloc( (n-1) * sizeof(float) );
    q_save = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    ab_r = (lapack_complex_float *)
        LAPACKE_malloc( (kd+1)*(n+2) * sizeof(lapack_complex_float) );
    q_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_ab( ldab*n, ab );
    init_d( n, d );
    init_e( (n-1), e );
    init_q( ldq*n, q );
    init_work( n, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldab*n; i++ ) {
        ab_save[i] = ab[i];
    }
    for( i = 0; i < n; i++ ) {
        d_save[i] = d[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_save[i] = e[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_save[i] = q[i];
    }

    /* Call the LAPACK routine */
    chbtrd_( &vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_chbtrd_work( LAPACK_COL_MAJOR, vect_i, uplo_i, n_i, kd_i,
                                  ab_i, ldab_i, d_i, e_i, q_i, ldq_i, work_i );

    failed = compare_chbtrd( ab, ab_i, d, d_i, e, e_i, q, q_i, info, info_i,
                             ldab, ldq, n, vect );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to chbtrd\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to chbtrd\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_chbtrd( LAPACK_COL_MAJOR, vect_i, uplo_i, n_i, kd_i, ab_i,
                             ldab_i, d_i, e_i, q_i, ldq_i );

    failed = compare_chbtrd( ab, ab_i, d, d_i, e, e_i, q, q_i, info, info_i,
                             ldab, ldq, n, vect );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to chbtrd\n" );
    } else {
        printf( "FAILED: column-major high-level interface to chbtrd\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, kd+1, n, ab_i, ldab, ab_r, n+2 );
    if( LAPACKE_lsame( vect, 'u' ) || LAPACKE_lsame( vect, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_chbtrd_work( LAPACK_ROW_MAJOR, vect_i, uplo_i, n_i, kd_i,
                                  ab_r, ldab_r, d_i, e_i, q_r, ldq_r, work_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, kd+1, n, ab_r, n+2, ab_i, ldab );
    if( LAPACKE_lsame( vect, 'u' ) || LAPACKE_lsame( vect, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_chbtrd( ab, ab_i, d, d_i, e, e_i, q, q_i, info, info_i,
                             ldab, ldq, n, vect );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to chbtrd\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to chbtrd\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldab*n; i++ ) {
        ab_i[i] = ab_save[i];
    }
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    for( i = 0; i < n; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, kd+1, n, ab_i, ldab, ab_r, n+2 );
    if( LAPACKE_lsame( vect, 'u' ) || LAPACKE_lsame( vect, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_chbtrd( LAPACK_ROW_MAJOR, vect_i, uplo_i, n_i, kd_i, ab_r,
                             ldab_r, d_i, e_i, q_r, ldq_r );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, kd+1, n, ab_r, n+2, ab_i, ldab );
    if( LAPACKE_lsame( vect, 'u' ) || LAPACKE_lsame( vect, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_chbtrd( ab, ab_i, d, d_i, e, e_i, q, q_i, info, info_i,
                             ldab, ldq, n, vect );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to chbtrd\n" );
    } else {
        printf( "FAILED: row-major high-level interface to chbtrd\n" );
    }

    /* Release memory */
    if( ab != NULL ) {
        LAPACKE_free( ab );
    }
    if( ab_i != NULL ) {
        LAPACKE_free( ab_i );
    }
    if( ab_r != NULL ) {
        LAPACKE_free( ab_r );
    }
    if( ab_save != NULL ) {
        LAPACKE_free( ab_save );
    }
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
    if( q != NULL ) {
        LAPACKE_free( q );
    }
    if( q_i != NULL ) {
        LAPACKE_free( q_i );
    }
    if( q_r != NULL ) {
        LAPACKE_free( q_r );
    }
    if( q_save != NULL ) {
        LAPACKE_free( q_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: chbtrd scalar parameters initialization */
static void init_scalars_chbtrd( char *vect, char *uplo, lapack_int *n,
                                 lapack_int *kd, lapack_int *ldab,
                                 lapack_int *ldq )
{
    *vect = 'V';
    *uplo = 'L';
    *n = 4;
    *kd = 2;
    *ldab = 9;
    *ldq = 8;

    return;
}

/* Auxiliary functions: chbtrd array parameters initialization */
static void init_ab( lapack_int size, lapack_complex_float *ab ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ab[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    ab[0] = lapack_make_complex_float( -3.130000114e+000, 0.000000000e+000 );
    ab[9] = lapack_make_complex_float( -1.909999967e+000, 0.000000000e+000 );
    ab[18] = lapack_make_complex_float( -2.869999886e+000, 0.000000000e+000 );
    ab[27] = lapack_make_complex_float( 5.000000000e-001, 0.000000000e+000 );
    ab[1] = lapack_make_complex_float( 1.940000057e+000, 2.099999905e+000 );
    ab[10] = lapack_make_complex_float( -8.199999928e-001, 8.899999857e-001 );
    ab[19] = lapack_make_complex_float( -2.099999905e+000, 1.599999964e-001 );
    ab[28] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[2] = lapack_make_complex_float( -3.400000095e+000, -2.500000000e-001 );
    ab[11] = lapack_make_complex_float( -6.700000167e-001, -3.400000036e-001 );
    ab[20] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    ab[29] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_d( lapack_int size, float *d ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        d[i] = 0;
    }
}
static void init_e( lapack_int size, float *e ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        e[i] = 0;
    }
}
static void init_q( lapack_int size, lapack_complex_float *q ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        q[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    q[0] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[8] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[16] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[24] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[9] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[17] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[25] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[10] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[18] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[26] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[19] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[27] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}

/* Auxiliary function: C interface to chbtrd results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_chbtrd( lapack_complex_float *ab, lapack_complex_float *ab_i,
                           float *d, float *d_i, float *e, float *e_i,
                           lapack_complex_float *q, lapack_complex_float *q_i,
                           lapack_int info, lapack_int info_i, lapack_int ldab,
                           lapack_int ldq, lapack_int n, char vect )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldab*n; i++ ) {
        failed += compare_complex_floats(ab[i],ab_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_floats(d[i],d_i[i]);
    }
    for( i = 0; i < (n-1); i++ ) {
        failed += compare_floats(e[i],e_i[i]);
    }
    if( LAPACKE_lsame( vect, 'u' ) || LAPACKE_lsame( vect, 'v' ) ) {
        for( i = 0; i < ldq*n; i++ ) {
            failed += compare_complex_floats(q[i],q_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
