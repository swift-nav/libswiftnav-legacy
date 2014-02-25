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
* dstein_1 is the test program for the C interface to LAPACK
* routine dstein
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

static void init_scalars_dstein( lapack_int *n, lapack_int *m,
                                 lapack_int *ldz );
static void init_d( lapack_int size, double *d );
static void init_e( lapack_int size, double *e );
static void init_w( lapack_int size, double *w );
static void init_iblock( lapack_int size, lapack_int *iblock );
static void init_isplit( lapack_int size, lapack_int *isplit );
static void init_z( lapack_int size, double *z );
static void init_work( lapack_int size, double *work );
static void init_iwork( lapack_int size, lapack_int *iwork );
static void init_ifailv( lapack_int size, lapack_int *ifailv );
static int compare_dstein( double *z, double *z_i, lapack_int *ifailv,
                           lapack_int *ifailv_i, lapack_int info,
                           lapack_int info_i, lapack_int ldz, lapack_int m );

int main(void)
{
    /* Local scalars */
    lapack_int n, n_i;
    lapack_int m, m_i;
    lapack_int ldz, ldz_i;
    lapack_int ldz_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *d = NULL, *d_i = NULL;
    double *e = NULL, *e_i = NULL;
    double *w = NULL, *w_i = NULL;
    lapack_int *iblock = NULL, *iblock_i = NULL;
    lapack_int *isplit = NULL, *isplit_i = NULL;
    double *z = NULL, *z_i = NULL;
    double *work = NULL, *work_i = NULL;
    lapack_int *iwork = NULL, *iwork_i = NULL;
    lapack_int *ifailv = NULL, *ifailv_i = NULL;
    double *z_save = NULL;
    lapack_int *ifailv_save = NULL;
    double *z_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dstein( &n, &m, &ldz );
    ldz_r = m+2;
    n_i = n;
    m_i = m;
    ldz_i = ldz;

    /* Allocate memory for the LAPACK routine arrays */
    d = (double *)LAPACKE_malloc( n * sizeof(double) );
    e = (double *)LAPACKE_malloc( n * sizeof(double) );
    w = (double *)LAPACKE_malloc( n * sizeof(double) );
    iblock = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    isplit = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    z = (double *)LAPACKE_malloc( ldz*m * sizeof(double) );
    work = (double *)LAPACKE_malloc( 5*n * sizeof(double) );
    iwork = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    ifailv = (lapack_int *)LAPACKE_malloc( m * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    d_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    e_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    w_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    iblock_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    isplit_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    z_i = (double *)LAPACKE_malloc( ldz*m * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( 5*n * sizeof(double) );
    iwork_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    ifailv_i = (lapack_int *)LAPACKE_malloc( m * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    z_save = (double *)LAPACKE_malloc( ldz*m * sizeof(double) );
    ifailv_save = (lapack_int *)LAPACKE_malloc( m * sizeof(lapack_int) );

    /* Allocate memory for the row-major arrays */
    z_r = (double *)LAPACKE_malloc( n*(m+2) * sizeof(double) );

    /* Initialize input arrays */
    init_d( n, d );
    init_e( n, e );
    init_w( n, w );
    init_iblock( n, iblock );
    init_isplit( n, isplit );
    init_z( ldz*m, z );
    init_work( 5*n, work );
    init_iwork( n, iwork );
    init_ifailv( m, ifailv );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldz*m; i++ ) {
        z_save[i] = z[i];
    }
    for( i = 0; i < m; i++ ) {
        ifailv_save[i] = ifailv[i];
    }

    /* Call the LAPACK routine */
    dstein_( &n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifailv,
             &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d[i];
    }
    for( i = 0; i < n; i++ ) {
        e_i[i] = e[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w[i];
    }
    for( i = 0; i < n; i++ ) {
        iblock_i[i] = iblock[i];
    }
    for( i = 0; i < n; i++ ) {
        isplit_i[i] = isplit[i];
    }
    for( i = 0; i < ldz*m; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < 5*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    for( i = 0; i < m; i++ ) {
        ifailv_i[i] = ifailv_save[i];
    }
    info_i = LAPACKE_dstein_work( LAPACK_COL_MAJOR, n_i, d_i, e_i, m_i, w_i,
                                  iblock_i, isplit_i, z_i, ldz_i, work_i,
                                  iwork_i, ifailv_i );

    failed = compare_dstein( z, z_i, ifailv, ifailv_i, info, info_i, ldz, m );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dstein\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dstein\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d[i];
    }
    for( i = 0; i < n; i++ ) {
        e_i[i] = e[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w[i];
    }
    for( i = 0; i < n; i++ ) {
        iblock_i[i] = iblock[i];
    }
    for( i = 0; i < n; i++ ) {
        isplit_i[i] = isplit[i];
    }
    for( i = 0; i < ldz*m; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < 5*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    for( i = 0; i < m; i++ ) {
        ifailv_i[i] = ifailv_save[i];
    }
    info_i = LAPACKE_dstein( LAPACK_COL_MAJOR, n_i, d_i, e_i, m_i, w_i,
                             iblock_i, isplit_i, z_i, ldz_i, ifailv_i );

    failed = compare_dstein( z, z_i, ifailv, ifailv_i, info, info_i, ldz, m );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dstein\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dstein\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d[i];
    }
    for( i = 0; i < n; i++ ) {
        e_i[i] = e[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w[i];
    }
    for( i = 0; i < n; i++ ) {
        iblock_i[i] = iblock[i];
    }
    for( i = 0; i < n; i++ ) {
        isplit_i[i] = isplit[i];
    }
    for( i = 0; i < ldz*m; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < 5*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    for( i = 0; i < m; i++ ) {
        ifailv_i[i] = ifailv_save[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, m, z_i, ldz, z_r, m+2 );
    info_i = LAPACKE_dstein_work( LAPACK_ROW_MAJOR, n_i, d_i, e_i, m_i, w_i,
                                  iblock_i, isplit_i, z_r, ldz_r, work_i,
                                  iwork_i, ifailv_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, m, z_r, m+2, z_i, ldz );

    failed = compare_dstein( z, z_i, ifailv, ifailv_i, info, info_i, ldz, m );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dstein\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dstein\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d[i];
    }
    for( i = 0; i < n; i++ ) {
        e_i[i] = e[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w[i];
    }
    for( i = 0; i < n; i++ ) {
        iblock_i[i] = iblock[i];
    }
    for( i = 0; i < n; i++ ) {
        isplit_i[i] = isplit[i];
    }
    for( i = 0; i < ldz*m; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < 5*n; i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < n; i++ ) {
        iwork_i[i] = iwork[i];
    }
    for( i = 0; i < m; i++ ) {
        ifailv_i[i] = ifailv_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, m, z_i, ldz, z_r, m+2 );
    info_i = LAPACKE_dstein( LAPACK_ROW_MAJOR, n_i, d_i, e_i, m_i, w_i,
                             iblock_i, isplit_i, z_r, ldz_r, ifailv_i );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, m, z_r, m+2, z_i, ldz );

    failed = compare_dstein( z, z_i, ifailv, ifailv_i, info, info_i, ldz, m );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dstein\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dstein\n" );
    }

    /* Release memory */
    if( d != NULL ) {
        LAPACKE_free( d );
    }
    if( d_i != NULL ) {
        LAPACKE_free( d_i );
    }
    if( e != NULL ) {
        LAPACKE_free( e );
    }
    if( e_i != NULL ) {
        LAPACKE_free( e_i );
    }
    if( w != NULL ) {
        LAPACKE_free( w );
    }
    if( w_i != NULL ) {
        LAPACKE_free( w_i );
    }
    if( iblock != NULL ) {
        LAPACKE_free( iblock );
    }
    if( iblock_i != NULL ) {
        LAPACKE_free( iblock_i );
    }
    if( isplit != NULL ) {
        LAPACKE_free( isplit );
    }
    if( isplit_i != NULL ) {
        LAPACKE_free( isplit_i );
    }
    if( z != NULL ) {
        LAPACKE_free( z );
    }
    if( z_i != NULL ) {
        LAPACKE_free( z_i );
    }
    if( z_r != NULL ) {
        LAPACKE_free( z_r );
    }
    if( z_save != NULL ) {
        LAPACKE_free( z_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }
    if( iwork != NULL ) {
        LAPACKE_free( iwork );
    }
    if( iwork_i != NULL ) {
        LAPACKE_free( iwork_i );
    }
    if( ifailv != NULL ) {
        LAPACKE_free( ifailv );
    }
    if( ifailv_i != NULL ) {
        LAPACKE_free( ifailv_i );
    }
    if( ifailv_save != NULL ) {
        LAPACKE_free( ifailv_save );
    }

    return 0;
}

/* Auxiliary function: dstein scalar parameters initialization */
static void init_scalars_dstein( lapack_int *n, lapack_int *m, lapack_int *ldz )
{
    *n = 4;
    *m = 2;
    *ldz = 8;

    return;
}

/* Auxiliary functions: dstein array parameters initialization */
static void init_d( lapack_int size, double *d ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        d[i] = 0;
    }
    d[0] = 2.06999999999999980e+000;
    d[1] = 1.47409370819755310e+000;
    d[2] = -6.49159507545784330e-001;
    d[3] = -1.69493420065176800e+000;
}
static void init_e( lapack_int size, double *e ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        e[i] = 0;
    }
    e[0] = -5.82575317019181590e+000;
    e[1] = 2.62404517879558740e+000;
    e[2] = 9.16272756321918620e-001;
    e[3] = 0.00000000000000000e+000;
}
static void init_w( lapack_int size, double *w ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        w[i] = 0;
    }
    w[0] = -5.00335314549203240e+000;
    w[1] = -1.99870457480866890e+000;
    w[2] = 0.00000000000000000e+000;
    w[3] = 0.00000000000000000e+000;
}
static void init_iblock( lapack_int size, lapack_int *iblock ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        iblock[i] = 0;
    }
    iblock[0] = 1;
    iblock[1] = 1;
    iblock[2] = 0;
    iblock[3] = 0;
}
static void init_isplit( lapack_int size, lapack_int *isplit ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        isplit[i] = 0;
    }
    isplit[0] = 4;
    isplit[1] = 0;
    isplit[2] = 0;
    isplit[3] = 0;
}
static void init_z( lapack_int size, double *z ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        z[i] = 0;
    }
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}
static void init_iwork( lapack_int size, lapack_int *iwork ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        iwork[i] = 0;
    }
}
static void init_ifailv( lapack_int size, lapack_int *ifailv ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ifailv[i] = 0;
    }
}

/* Auxiliary function: C interface to dstein results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dstein( double *z, double *z_i, lapack_int *ifailv,
                           lapack_int *ifailv_i, lapack_int info,
                           lapack_int info_i, lapack_int ldz, lapack_int m )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldz*m; i++ ) {
        failed += compare_doubles(z[i],z_i[i]);
    }
    for( i = 0; i < m; i++ ) {
        failed += (ifailv[i] == ifailv_i[i]) ? 0 : 1;
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
