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
* dsteqr_1 is the test program for the C interface to LAPACK
* routine dsteqr
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

static void init_scalars_dsteqr( char *compz, lapack_int *n, lapack_int *ldz );
static void init_d( lapack_int size, double *d );
static void init_e( lapack_int size, double *e );
static void init_z( lapack_int size, double *z );
static void init_work( lapack_int size, double *work );
static int compare_dsteqr( double *d, double *d_i, double *e, double *e_i,
                           double *z, double *z_i, lapack_int info,
                           lapack_int info_i, char compz, lapack_int ldz,
                           lapack_int n );

int main(void)
{
    /* Local scalars */
    char compz, compz_i;
    lapack_int n, n_i;
    lapack_int ldz, ldz_i;
    lapack_int ldz_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *d = NULL, *d_i = NULL;
    double *e = NULL, *e_i = NULL;
    double *z = NULL, *z_i = NULL;
    double *work = NULL, *work_i = NULL;
    double *d_save = NULL;
    double *e_save = NULL;
    double *z_save = NULL;
    double *z_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dsteqr( &compz, &n, &ldz );
    ldz_r = n+2;
    compz_i = compz;
    n_i = n;
    ldz_i = ldz;

    /* Allocate memory for the LAPACK routine arrays */
    d = (double *)LAPACKE_malloc( n * sizeof(double) );
    e = (double *)LAPACKE_malloc( (n-1) * sizeof(double) );
    z = (double *)LAPACKE_malloc( ldz*n * sizeof(double) );
    work = (double *)LAPACKE_malloc( ((MAX(1,2*n-2))) * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    d_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    e_i = (double *)LAPACKE_malloc( (n-1) * sizeof(double) );
    z_i = (double *)LAPACKE_malloc( ldz*n * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( ((MAX(1,2*n-2))) * sizeof(double) );

    /* Allocate memory for the backup arrays */
    d_save = (double *)LAPACKE_malloc( n * sizeof(double) );
    e_save = (double *)LAPACKE_malloc( (n-1) * sizeof(double) );
    z_save = (double *)LAPACKE_malloc( ldz*n * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    z_r = (double *)LAPACKE_malloc( n*(n+2) * sizeof(double) );

    /* Initialize input arrays */
    init_d( n, d );
    init_e( (n-1), e );
    init_z( ldz*n, z );
    init_work( (MAX(1,2*n-2)), work );

    /* Backup the ouptut arrays */
    for( i = 0; i < n; i++ ) {
        d_save[i] = d[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_save[i] = e[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_save[i] = z[i];
    }

    /* Call the LAPACK routine */
    dsteqr_( &compz, &n, d, e, z, &ldz, work, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < (MAX(1,2*n-2)); i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dsteqr_work( LAPACK_COL_MAJOR, compz_i, n_i, d_i, e_i, z_i,
                                  ldz_i, work_i );

    failed = compare_dsteqr( d, d_i, e, e_i, z, z_i, info, info_i, compz, ldz,
                             n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dsteqr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dsteqr\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < (MAX(1,2*n-2)); i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_dsteqr( LAPACK_COL_MAJOR, compz_i, n_i, d_i, e_i, z_i,
                             ldz_i );

    failed = compare_dsteqr( d, d_i, e, e_i, z, z_i, info, info_i, compz, ldz,
                             n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dsteqr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dsteqr\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < (MAX(1,2*n-2)); i++ ) {
        work_i[i] = work[i];
    }

    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, z_i, ldz, z_r, n+2 );
    }
    info_i = LAPACKE_dsteqr_work( LAPACK_ROW_MAJOR, compz_i, n_i, d_i, e_i, z_r,
                                  ldz_r, work_i );

    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, n, z_r, n+2, z_i, ldz );
    }

    failed = compare_dsteqr( d, d_i, e, e_i, z, z_i, info, info_i, compz, ldz,
                             n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dsteqr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dsteqr\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        d_i[i] = d_save[i];
    }
    for( i = 0; i < (n-1); i++ ) {
        e_i[i] = e_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < (MAX(1,2*n-2)); i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, z_i, ldz, z_r, n+2 );
    }
    info_i = LAPACKE_dsteqr( LAPACK_ROW_MAJOR, compz_i, n_i, d_i, e_i, z_r,
                             ldz_r );

    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, n, z_r, n+2, z_i, ldz );
    }

    failed = compare_dsteqr( d, d_i, e, e_i, z, z_i, info, info_i, compz, ldz,
                             n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dsteqr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dsteqr\n" );
    }

    /* Release memory */
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

    return 0;
}

/* Auxiliary function: dsteqr scalar parameters initialization */
static void init_scalars_dsteqr( char *compz, lapack_int *n, lapack_int *ldz )
{
    *compz = 'V';
    *n = 4;
    *ldz = 8;

    return;
}

/* Auxiliary functions: dsteqr array parameters initialization */
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
}
static void init_z( lapack_int size, double *z ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        z[i] = 0;
    }
    z[0] = 1.00000000000000000e+000;  /* z[0,0] */
    z[8] = 0.00000000000000000e+000;  /* z[0,1] */
    z[16] = 0.00000000000000000e+000;  /* z[0,2] */
    z[24] = 0.00000000000000000e+000;  /* z[0,3] */
    z[1] = 0.00000000000000000e+000;  /* z[1,0] */
    z[9] = -6.64291789738249210e-001;  /* z[1,1] */
    z[17] = -4.00376638335342840e-002;  /* z[1,2] */
    z[25] = 7.46400297133587130e-001;  /* z[1,3] */
    z[2] = 0.00000000000000000e+000;  /* z[2,0] */
    z[10] = -7.20936826072518520e-001;  /* z[2,1] */
    z[18] = -2.29390813125724850e-001;  /* z[2,2] */
    z[26] = -6.53934207444450740e-001;  /* z[2,3] */
    z[3] = 0.00000000000000000e+000;  /* z[3,0] */
    z[11] = 1.97399369043665780e-001;  /* z[3,1] */
    z[19] = -9.72510586229410450e-001;  /* z[3,2] */
    z[27] = 1.23517807513293280e-001;  /* z[3,3] */
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to dsteqr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dsteqr( double *d, double *d_i, double *e, double *e_i,
                           double *z, double *z_i, lapack_int info,
                           lapack_int info_i, char compz, lapack_int ldz,
                           lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < n; i++ ) {
        failed += compare_doubles(d[i],d_i[i]);
    }
    for( i = 0; i < (n-1); i++ ) {
        failed += compare_doubles(e[i],e_i[i]);
    }
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        for( i = 0; i < ldz*n; i++ ) {
            failed += compare_doubles(z[i],z_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
