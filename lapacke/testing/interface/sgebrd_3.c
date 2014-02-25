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
* sgebrd_3 is the test program for the C interface to LAPACK
* routine sgebrd
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

static void init_scalars_sgebrd( lapack_int *m, lapack_int *n, lapack_int *lda,
                                 lapack_int *lwork );
static void init_a( lapack_int size, float *a );
static void init_d( lapack_int size, float *d );
static void init_e( lapack_int size, float *e );
static void init_tauq( lapack_int size, float *tauq );
static void init_taup( lapack_int size, float *taup );
static void init_work( lapack_int size, float *work );
static int compare_sgebrd( float *a, float *a_i, float *d, float *d_i, float *e,
                           float *e_i, float *tauq, float *tauq_i, float *taup,
                           float *taup_i, lapack_int info, lapack_int info_i,
                           lapack_int lda, lapack_int m, lapack_int n );

int main(void)
{
    /* Local scalars */
    lapack_int m, m_i;
    lapack_int n, n_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    lapack_int lwork, lwork_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    float *a = NULL, *a_i = NULL;
    float *d = NULL, *d_i = NULL;
    float *e = NULL, *e_i = NULL;
    float *tauq = NULL, *tauq_i = NULL;
    float *taup = NULL, *taup_i = NULL;
    float *work = NULL, *work_i = NULL;
    float *a_save = NULL;
    float *d_save = NULL;
    float *e_save = NULL;
    float *tauq_save = NULL;
    float *taup_save = NULL;
    float *a_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_sgebrd( &m, &n, &lda, &lwork );
    lda_r = n+2;
    m_i = m;
    n_i = n;
    lda_i = lda;
    lwork_i = lwork;

    /* Allocate memory for the LAPACK routine arrays */
    a = (float *)LAPACKE_malloc( lda*n * sizeof(float) );
    d = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );
    e = (float *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(float) );
    tauq = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );
    taup = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );
    work = (float *)LAPACKE_malloc( lwork * sizeof(float) );

    /* Allocate memory for the C interface function arrays */
    a_i = (float *)LAPACKE_malloc( lda*n * sizeof(float) );
    d_i = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );
    e_i = (float *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(float) );
    tauq_i = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );
    taup_i = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );
    work_i = (float *)LAPACKE_malloc( lwork * sizeof(float) );

    /* Allocate memory for the backup arrays */
    a_save = (float *)LAPACKE_malloc( lda*n * sizeof(float) );
    d_save = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );
    e_save = (float *)LAPACKE_malloc( ((MIN(m,n)-1)) * sizeof(float) );
    tauq_save = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );
    taup_save = (float *)LAPACKE_malloc( MIN(m,n) * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    a_r = (float *)LAPACKE_malloc( m*(n+2) * sizeof(float) );

    /* Initialize input arrays */
    init_a( lda*n, a );
    init_d( (MIN(m,n)), d );
    init_e( (MIN(m,n)-1), e );
    init_tauq( (MIN(m,n)), tauq );
    init_taup( (MIN(m,n)), taup );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < lda*n; i++ ) {
        a_save[i] = a[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_save[i] = d[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_save[i] = e[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tauq_save[i] = tauq[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        taup_save[i] = taup[i];
    }

    /* Call the LAPACK routine */
    sgebrd_( &m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tauq_i[i] = tauq_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        taup_i[i] = taup_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_sgebrd_work( LAPACK_COL_MAJOR, m_i, n_i, a_i, lda_i, d_i,
                                  e_i, tauq_i, taup_i, work_i, lwork_i );

    failed = compare_sgebrd( a, a_i, d, d_i, e, e_i, tauq, tauq_i, taup, taup_i,
                             info, info_i, lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to sgebrd\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to sgebrd\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tauq_i[i] = tauq_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        taup_i[i] = taup_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_sgebrd( LAPACK_COL_MAJOR, m_i, n_i, a_i, lda_i, d_i, e_i,
                             tauq_i, taup_i );

    failed = compare_sgebrd( a, a_i, d, d_i, e, e_i, tauq, tauq_i, taup, taup_i,
                             info, info_i, lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to sgebrd\n" );
    } else {
        printf( "FAILED: column-major high-level interface to sgebrd\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tauq_i[i] = tauq_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        taup_i[i] = taup_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_sgebrd_work( LAPACK_ROW_MAJOR, m_i, n_i, a_r, lda_r, d_i,
                                  e_i, tauq_i, taup_i, work_i, lwork_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, m, n, a_r, n+2, a_i, lda );

    failed = compare_sgebrd( a, a_i, d, d_i, e, e_i, tauq, tauq_i, taup, taup_i,
                             info, info_i, lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to sgebrd\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to sgebrd\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*n; i++ ) {
        a_i[i] = a_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        tauq_i[i] = tauq_save[i];
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        taup_i[i] = taup_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, m, n, a_i, lda, a_r, n+2 );
    info_i = LAPACKE_sgebrd( LAPACK_ROW_MAJOR, m_i, n_i, a_r, lda_r, d_i, e_i,
                             tauq_i, taup_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, m, n, a_r, n+2, a_i, lda );

    failed = compare_sgebrd( a, a_i, d, d_i, e, e_i, tauq, tauq_i, taup, taup_i,
                             info, info_i, lda, m, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to sgebrd\n" );
    } else {
        printf( "FAILED: row-major high-level interface to sgebrd\n" );
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
    if( tauq != NULL ) {
        LAPACKE_free( tauq );
    }
    if( tauq_i != NULL ) {
        LAPACKE_free( tauq_i );
    }
    if( tauq_save != NULL ) {
        LAPACKE_free( tauq_save );
    }
    if( taup != NULL ) {
        LAPACKE_free( taup );
    }
    if( taup_i != NULL ) {
        LAPACKE_free( taup_i );
    }
    if( taup_save != NULL ) {
        LAPACKE_free( taup_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }

    return 0;
}

/* Auxiliary function: sgebrd scalar parameters initialization */
static void init_scalars_sgebrd( lapack_int *m, lapack_int *n, lapack_int *lda,
                                 lapack_int *lwork )
{
    *m = 4;
    *n = 6;
    *lda = 8;
    *lwork = 1024;

    return;
}

/* Auxiliary functions: sgebrd array parameters initialization */
static void init_a( lapack_int size, float *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = 0;
    }
    a[0] = -5.420000076e+000;  /* a[0,0] */
    a[8] = 3.279999971e+000;  /* a[0,1] */
    a[16] = -3.680000067e+000;  /* a[0,2] */
    a[24] = 2.700000107e-001;  /* a[0,3] */
    a[32] = 2.059999943e+000;  /* a[0,4] */
    a[40] = 4.600000083e-001;  /* a[0,5] */
    a[1] = -1.649999976e+000;  /* a[1,0] */
    a[9] = -3.400000095e+000;  /* a[1,1] */
    a[17] = -3.200000048e+000;  /* a[1,2] */
    a[25] = -1.029999971e+000;  /* a[1,3] */
    a[33] = -4.059999943e+000;  /* a[1,4] */
    a[41] = -9.999999776e-003;  /* a[1,5] */
    a[2] = -3.700000048e-001;  /* a[2,0] */
    a[10] = 2.349999905e+000;  /* a[2,1] */
    a[18] = 1.899999976e+000;  /* a[2,2] */
    a[26] = 4.309999943e+000;  /* a[2,3] */
    a[34] = -1.759999990e+000;  /* a[2,4] */
    a[42] = 1.129999995e+000;  /* a[2,5] */
    a[3] = -3.150000095e+000;  /* a[3,0] */
    a[11] = -1.099999994e-001;  /* a[3,1] */
    a[19] = 1.990000010e+000;  /* a[3,2] */
    a[27] = -2.700000048e+000;  /* a[3,3] */
    a[35] = 2.599999905e-001;  /* a[3,4] */
    a[43] = 4.500000000e+000;  /* a[3,5] */
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
static void init_tauq( lapack_int size, float *tauq ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tauq[i] = 0;
    }
}
static void init_taup( lapack_int size, float *taup ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        taup[i] = 0;
    }
}
static void init_work( lapack_int size, float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to sgebrd results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_sgebrd( float *a, float *a_i, float *d, float *d_i, float *e,
                           float *e_i, float *tauq, float *tauq_i, float *taup,
                           float *taup_i, lapack_int info, lapack_int info_i,
                           lapack_int lda, lapack_int m, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < lda*n; i++ ) {
        failed += compare_floats(a[i],a_i[i]);
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        failed += compare_floats(d[i],d_i[i]);
    }
    for( i = 0; i < (MIN(m,n)-1); i++ ) {
        failed += compare_floats(e[i],e_i[i]);
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        failed += compare_floats(tauq[i],tauq_i[i]);
    }
    for( i = 0; i < (MIN(m,n)); i++ ) {
        failed += compare_floats(taup[i],taup_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
