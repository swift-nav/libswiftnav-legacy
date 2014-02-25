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
* dormbr_1 is the test program for the C interface to LAPACK
* routine dormbr
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

static void init_scalars_dormbr( char *vect, char *side, char *trans,
                                 lapack_int *m, lapack_int *n, lapack_int *k,
                                 lapack_int *lda, lapack_int *ldc,
                                 lapack_int *lwork );
static void init_a( lapack_int size, double *a );
static void init_tau( lapack_int size, double *tau );
static void init_c( lapack_int size, double *c );
static void init_work( lapack_int size, double *work );
static int compare_dormbr( double *c, double *c_i, lapack_int info,
                           lapack_int info_i, lapack_int ldc, lapack_int n );

int main(void)
{
    /* Local scalars */
    char vect, vect_i;
    char side, side_i;
    char trans, trans_i;
    lapack_int m, m_i;
    lapack_int n, n_i;
    lapack_int k, k_i;
    lapack_int lda, lda_i;
    lapack_int lda_r;
    lapack_int ldc, ldc_i;
    lapack_int ldc_r;
    lapack_int lwork, lwork_i;
    lapack_int info, info_i;
    /* Declare scalars */
    lapack_int nq;
    lapack_int r;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *a = NULL, *a_i = NULL;
    double *tau = NULL, *tau_i = NULL;
    double *c = NULL, *c_i = NULL;
    double *work = NULL, *work_i = NULL;
    double *c_save = NULL;
    double *a_r = NULL;
    double *c_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dormbr( &vect, &side, &trans, &m, &n, &k, &lda, &ldc, &lwork );
    nq = LAPACKE_lsame( side, 'l' ) ? m : n;
    r = LAPACKE_lsame( vect, 'q' ) ? nq : MIN(nq,k);
    lda_r = MIN(nq,k)+2;
    ldc_r = n+2;
    vect_i = vect;
    side_i = side;
    trans_i = trans;
    m_i = m;
    n_i = n;
    k_i = k;
    lda_i = lda;
    ldc_i = ldc;
    lwork_i = lwork;

    /* Allocate memory for the LAPACK routine arrays */
    a = (double *)LAPACKE_malloc( (lda*(MIN(nq,k))) * sizeof(double) );
    tau = (double *)LAPACKE_malloc( MIN(nq,k) * sizeof(double) );
    c = (double *)LAPACKE_malloc( ldc*n * sizeof(double) );
    work = (double *)LAPACKE_malloc( lwork * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    a_i = (double *)LAPACKE_malloc( (lda*(MIN(nq,k))) * sizeof(double) );
    tau_i = (double *)LAPACKE_malloc( MIN(nq,k) * sizeof(double) );
    c_i = (double *)LAPACKE_malloc( ldc*n * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( lwork * sizeof(double) );

    /* Allocate memory for the backup arrays */
    c_save = (double *)LAPACKE_malloc( ldc*n * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    a_r = (double *)LAPACKE_malloc( (r*(MIN(nq,k)+2)) * sizeof(double) );
    c_r = (double *)LAPACKE_malloc( m*(n+2) * sizeof(double) );

    /* Initialize input arrays */
    init_a( lda*(MIN(nq,k)), a );
    init_tau( (MIN(nq,k)), tau );
    init_c( ldc*n, c );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldc*n; i++ ) {
        c_save[i] = c[i];
    }

    /* Call the LAPACK routine */
    dormbr_( &vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work,
             &lwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*(MIN(nq,k)); i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < (MIN(nq,k)); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dormbr_work( LAPACK_COL_MAJOR, vect_i, side_i, trans_i,
                                  m_i, n_i, k_i, a_i, lda_i, tau_i, c_i, ldc_i,
                                  work_i, lwork_i );

    failed = compare_dormbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dormbr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dormbr\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*(MIN(nq,k)); i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < (MIN(nq,k)); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dormbr( LAPACK_COL_MAJOR, vect_i, side_i, trans_i, m_i,
                             n_i, k_i, a_i, lda_i, tau_i, c_i, ldc_i );

    failed = compare_dormbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dormbr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dormbr\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*(MIN(nq,k)); i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < (MIN(nq,k)); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, r, MIN(nq, k ), a_i, lda, a_r, MIN(nq,
                       k)+2);
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_dormbr_work( LAPACK_ROW_MAJOR, vect_i, side_i, trans_i,
                                  m_i, n_i, k_i, a_r, lda_r, tau_i, c_r, ldc_r,
                                  work_i, lwork_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_dormbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dormbr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dormbr\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < lda*(MIN(nq,k)); i++ ) {
        a_i[i] = a[i];
    }
    for( i = 0; i < (MIN(nq,k)); i++ ) {
        tau_i[i] = tau[i];
    }
    for( i = 0; i < ldc*n; i++ ) {
        c_i[i] = c_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, r, MIN(nq, k ), a_i, lda, a_r, MIN(nq,
                       k)+2);
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, m, n, c_i, ldc, c_r, n+2 );
    info_i = LAPACKE_dormbr( LAPACK_ROW_MAJOR, vect_i, side_i, trans_i, m_i,
                             n_i, k_i, a_r, lda_r, tau_i, c_r, ldc_r );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, m, n, c_r, n+2, c_i, ldc );

    failed = compare_dormbr( c, c_i, info, info_i, ldc, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dormbr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dormbr\n" );
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
    if( tau != NULL ) {
        LAPACKE_free( tau );
    }
    if( tau_i != NULL ) {
        LAPACKE_free( tau_i );
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

/* Auxiliary function: dormbr scalar parameters initialization */
static void init_scalars_dormbr( char *vect, char *side, char *trans,
                                 lapack_int *m, lapack_int *n, lapack_int *k,
                                 lapack_int *lda, lapack_int *ldc,
                                 lapack_int *lwork )
{
    *vect = 'Q';
    *side = 'R';
    *trans = 'N';
    *m = 6;
    *n = 4;
    *k = 4;
    *lda = 8;
    *ldc = 8;
    *lwork = 1024;

    return;
}

/* Auxiliary functions: dormbr array parameters initialization */
static void init_a( lapack_int size, double *a ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        a[i] = 0;
    }
    a[0] = 3.61767881382524030e+000;  /* a[0,0] */
    a[8] = 1.25871140471055650e+000;  /* a[0,1] */
    a[16] = -4.66788533970659230e-001;  /* a[0,2] */
    a[24] = -4.10950544929054030e-001;  /* a[0,3] */
    a[1] = 0.00000000000000000e+000;  /* a[1,0] */
    a[9] = -2.41605534042332600e+000;  /* a[1,1] */
    a[17] = -1.52615488994111680e+000;  /* a[1,2] */
    a[25] = -2.09455650564123440e-001;  /* a[1,3] */
    a[2] = 0.00000000000000000e+000;  /* a[2,0] */
    a[10] = 2.04091961983997000e-002;  /* a[2,1] */
    a[18] = 1.92131096370999940e+000;  /* a[2,2] */
    a[26] = 1.18946667701939000e+000;  /* a[2,3] */
    a[3] = 0.00000000000000000e+000;  /* a[3,0] */
    a[11] = -3.21626418690548640e-001;  /* a[3,1] */
    a[19] = 1.34228166612422040e-001;  /* a[3,2] */
    a[27] = -1.42650154173085970e+000;  /* a[3,3] */
}
static void init_tau( lapack_int size, double *tau ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        tau[i] = 0;
    }
    tau[0] = 0.00000000000000000e+000;
    tau[1] = 1.81182381794177600e+000;
    tau[2] = 1.96460334717423550e+000;
    tau[3] = 0.00000000000000000e+000;
}
static void init_c( lapack_int size, double *c ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        c[i] = 0;
    }
    c[0] = -1.57559592582321220e-001;  /* c[0,0] */
    c[8] = 6.74381518209104080e-001;  /* c[0,1] */
    c[16] = -4.57149964215295970e-001;  /* c[0,2] */
    c[24] = 4.48851645119849540e-001;  /* c[0,3] */
    c[1] = -5.33491252076982340e-001;  /* c[1,0] */
    c[9] = -3.86108989331406360e-001;  /* c[1,1] */
    c[17] = 2.58252726077486450e-001;  /* c[1,2] */
    c[25] = 3.89817437867655880e-001;  /* c[1,3] */
    c[2] = 6.35766777086559380e-001;  /* c[2,0] */
    c[10] = -2.92823089959412750e-001;  /* c[2,1] */
    c[18] = 1.65384263845584910e-002;  /* c[2,2] */
    c[26] = 1.92952916776812260e-001;  /* c[2,3] */
    c[3] = -5.33491252076982340e-001;  /* c[3,0] */
    c[11] = -1.69154705537664260e-001;  /* c[3,1] */
    c[19] = -8.34273545141680510e-002;  /* c[3,2] */
    c[27] = -2.34987319363205350e-001;  /* c[3,3] */
    c[4] = 4.14630506795582170e-002;  /* c[4,0] */
    c[12] = -1.59302792034020260e-001;  /* c[4,1] */
    c[20] = 1.47478303793470600e-001;  /* c[4,2] */
    c[28] = 7.43640764641989470e-001;  /* c[4,3] */
    c[5] = -5.52840675727442940e-003;  /* c[5,0] */
    c[13] = -5.06352999681760840e-001;  /* c[5,1] */
    c[21] = -8.33868063049059920e-001;  /* c[5,2] */
    c[29] = 3.35128425531637520e-002;  /* c[5,3] */
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to dormbr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dormbr( double *c, double *c_i, lapack_int info,
                           lapack_int info_i, lapack_int ldc, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldc*n; i++ ) {
        failed += compare_doubles(c[i],c_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
