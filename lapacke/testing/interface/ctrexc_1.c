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
* ctrexc_1 is the test program for the C interface to LAPACK
* routine ctrexc
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

static void init_scalars_ctrexc( char *compq, lapack_int *n, lapack_int *ldt,
                                 lapack_int *ldq, lapack_int *ifst,
                                 lapack_int *ilst );
static void init_t( lapack_int size, lapack_complex_float *t );
static void init_q( lapack_int size, lapack_complex_float *q );
static int compare_ctrexc( lapack_complex_float *t, lapack_complex_float *t_i,
                           lapack_complex_float *q, lapack_complex_float *q_i,
                           lapack_int info, lapack_int info_i, char compq,
                           lapack_int ldq, lapack_int ldt, lapack_int n );

int main(void)
{
    /* Local scalars */
    char compq, compq_i;
    lapack_int n, n_i;
    lapack_int ldt, ldt_i;
    lapack_int ldt_r;
    lapack_int ldq, ldq_i;
    lapack_int ldq_r;
    lapack_int ifst, ifst_i;
    lapack_int ilst, ilst_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_complex_float *t = NULL, *t_i = NULL;
    lapack_complex_float *q = NULL, *q_i = NULL;
    lapack_complex_float *t_save = NULL;
    lapack_complex_float *q_save = NULL;
    lapack_complex_float *t_r = NULL;
    lapack_complex_float *q_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_ctrexc( &compq, &n, &ldt, &ldq, &ifst, &ilst );
    ldt_r = n+2;
    ldq_r = n+2;
    compq_i = compq;
    n_i = n;
    ldt_i = ldt;
    ldq_i = ldq;
    ifst_i = ifst;
    ilst_i = ilst;

    /* Allocate memory for the LAPACK routine arrays */
    t = (lapack_complex_float *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_float) );
    q = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    t_i = (lapack_complex_float *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_float) );
    q_i = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    t_save = (lapack_complex_float *)
        LAPACKE_malloc( ldt*n * sizeof(lapack_complex_float) );
    q_save = (lapack_complex_float *)
        LAPACKE_malloc( ldq*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    t_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );
    q_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );

    /* Initialize input arrays */
    init_t( ldt*n, t );
    init_q( ldq*n, q );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldt*n; i++ ) {
        t_save[i] = t[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_save[i] = q[i];
    }

    /* Call the LAPACK routine */
    ctrexc_( &compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    info_i = LAPACKE_ctrexc_work( LAPACK_COL_MAJOR, compq_i, n_i, t_i, ldt_i,
                                  q_i, ldq_i, ifst_i, ilst_i );

    failed = compare_ctrexc( t, t_i, q, q_i, info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to ctrexc\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to ctrexc\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }
    info_i = LAPACKE_ctrexc( LAPACK_COL_MAJOR, compq_i, n_i, t_i, ldt_i, q_i,
                             ldq_i, ifst_i, ilst_i );

    failed = compare_ctrexc( t, t_i, q, q_i, info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to ctrexc\n" );
    } else {
        printf( "FAILED: column-major high-level interface to ctrexc\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_ctrexc_work( LAPACK_ROW_MAJOR, compq_i, n_i, t_r, ldt_r,
                                  q_r, ldq_r, ifst_i, ilst_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, t_r, n+2, t_i, ldt );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_ctrexc( t, t_i, q, q_i, info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to ctrexc\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to ctrexc\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldt*n; i++ ) {
        t_i[i] = t_save[i];
    }
    for( i = 0; i < ldq*n; i++ ) {
        q_i[i] = q_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, t_i, ldt, t_r, n+2 );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, q_i, ldq, q_r, n+2 );
    }
    info_i = LAPACKE_ctrexc( LAPACK_ROW_MAJOR, compq_i, n_i, t_r, ldt_r, q_r,
                             ldq_r, ifst_i, ilst_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, t_r, n+2, t_i, ldt );
    if( LAPACKE_lsame( compq, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, q_r, n+2, q_i, ldq );
    }

    failed = compare_ctrexc( t, t_i, q, q_i, info, info_i, compq, ldq, ldt, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to ctrexc\n" );
    } else {
        printf( "FAILED: row-major high-level interface to ctrexc\n" );
    }

    /* Release memory */
    if( t != NULL ) {
        LAPACKE_free( t );
    }
    if( t_i != NULL ) {
        LAPACKE_free( t_i );
    }
    if( t_r != NULL ) {
        LAPACKE_free( t_r );
    }
    if( t_save != NULL ) {
        LAPACKE_free( t_save );
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

    return 0;
}

/* Auxiliary function: ctrexc scalar parameters initialization */
static void init_scalars_ctrexc( char *compq, lapack_int *n, lapack_int *ldt,
                                 lapack_int *ldq, lapack_int *ifst,
                                 lapack_int *ilst )
{
    *compq = 'N';
    *n = 4;
    *ldt = 8;
    *ldq = 1;
    *ifst = 1;
    *ilst = 4;

    return;
}

/* Auxiliary functions: ctrexc array parameters initialization */
static void init_t( lapack_int size, lapack_complex_float *t ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        t[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    t[0] = lapack_make_complex_float( -6.000000000e+000, -7.000000000e+000 );
    t[8] = lapack_make_complex_float( 3.600000143e-001, -3.600000143e-001 );
    t[16] = lapack_make_complex_float( -1.899999976e-001, 4.799999893e-001 );
    t[24] = lapack_make_complex_float( 8.799999952e-001, -2.500000000e-001 );
    t[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[9] = lapack_make_complex_float( -5.000000000e+000, 2.000000000e+000 );
    t[17] = lapack_make_complex_float( -2.999999933e-002, -7.200000286e-001 );
    t[25] = lapack_make_complex_float( -2.300000042e-001, 1.299999952e-001 );
    t[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[10] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[18] = lapack_make_complex_float( 8.000000000e+000, -1.000000000e+000 );
    t[26] = lapack_make_complex_float( 9.399999976e-001, 5.299999714e-001 );
    t[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[11] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[19] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    t[27] = lapack_make_complex_float( 3.000000000e+000, -4.000000000e+000 );
}
static void init_q( lapack_int size, lapack_complex_float *q ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        q[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    q[0] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    q[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
}

/* Auxiliary function: C interface to ctrexc results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_ctrexc( lapack_complex_float *t, lapack_complex_float *t_i,
                           lapack_complex_float *q, lapack_complex_float *q_i,
                           lapack_int info, lapack_int info_i, char compq,
                           lapack_int ldq, lapack_int ldt, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldt*n; i++ ) {
        failed += compare_complex_floats(t[i],t_i[i]);
    }
    if( LAPACKE_lsame( compq, 'v' ) ) {
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
