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
* shseqr_2 is the test program for the C interface to LAPACK
* routine shseqr
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

static void init_scalars_shseqr( char *job, char *compz, lapack_int *n,
                                 lapack_int *ilo, lapack_int *ihi,
                                 lapack_int *ldh, lapack_int *ldz,
                                 lapack_int *lwork );
static void init_h( lapack_int size, float *h );
static void init_wr( lapack_int size, float *wr );
static void init_wi( lapack_int size, float *wi );
static void init_z( lapack_int size, float *z );
static void init_work( lapack_int size, float *work );
static int compare_shseqr( float *h, float *h_i, float *wr, float *wr_i,
                           float *wi, float *wi_i, float *z, float *z_i,
                           lapack_int info, lapack_int info_i, char compz,
                           lapack_int ldh, lapack_int ldz, lapack_int n );

int main(void)
{
    /* Local scalars */
    char job, job_i;
    char compz, compz_i;
    lapack_int n, n_i;
    lapack_int ilo, ilo_i;
    lapack_int ihi, ihi_i;
    lapack_int ldh, ldh_i;
    lapack_int ldh_r;
    lapack_int ldz, ldz_i;
    lapack_int ldz_r;
    lapack_int lwork, lwork_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    float *h = NULL, *h_i = NULL;
    float *wr = NULL, *wr_i = NULL;
    float *wi = NULL, *wi_i = NULL;
    float *z = NULL, *z_i = NULL;
    float *work = NULL, *work_i = NULL;
    float *h_save = NULL;
    float *wr_save = NULL;
    float *wi_save = NULL;
    float *z_save = NULL;
    float *h_r = NULL;
    float *z_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_shseqr( &job, &compz, &n, &ilo, &ihi, &ldh, &ldz, &lwork );
    ldh_r = n+2;
    ldz_r = n+2;
    job_i = job;
    compz_i = compz;
    n_i = n;
    ilo_i = ilo;
    ihi_i = ihi;
    ldh_i = ldh;
    ldz_i = ldz;
    lwork_i = lwork;

    /* Allocate memory for the LAPACK routine arrays */
    h = (float *)LAPACKE_malloc( ldh*n * sizeof(float) );
    wr = (float *)LAPACKE_malloc( n * sizeof(float) );
    wi = (float *)LAPACKE_malloc( n * sizeof(float) );
    z = (float *)LAPACKE_malloc( ldz*n * sizeof(float) );
    work = (float *)LAPACKE_malloc( lwork * sizeof(float) );

    /* Allocate memory for the C interface function arrays */
    h_i = (float *)LAPACKE_malloc( ldh*n * sizeof(float) );
    wr_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    wi_i = (float *)LAPACKE_malloc( n * sizeof(float) );
    z_i = (float *)LAPACKE_malloc( ldz*n * sizeof(float) );
    work_i = (float *)LAPACKE_malloc( lwork * sizeof(float) );

    /* Allocate memory for the backup arrays */
    h_save = (float *)LAPACKE_malloc( ldh*n * sizeof(float) );
    wr_save = (float *)LAPACKE_malloc( n * sizeof(float) );
    wi_save = (float *)LAPACKE_malloc( n * sizeof(float) );
    z_save = (float *)LAPACKE_malloc( ldz*n * sizeof(float) );

    /* Allocate memory for the row-major arrays */
    h_r = (float *)LAPACKE_malloc( n*(n+2) * sizeof(float) );
    z_r = (float *)LAPACKE_malloc( n*(n+2) * sizeof(float) );

    /* Initialize input arrays */
    init_h( ldh*n, h );
    init_wr( n, wr );
    init_wi( n, wi );
    init_z( ldz*n, z );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldh*n; i++ ) {
        h_save[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_save[i] = wr[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_save[i] = wi[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_save[i] = z[i];
    }

    /* Call the LAPACK routine */
    shseqr_( &job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, work,
             &lwork, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_shseqr_work( LAPACK_COL_MAJOR, job_i, compz_i, n_i, ilo_i,
                                  ihi_i, h_i, ldh_i, wr_i, wi_i, z_i, ldz_i,
                                  work_i, lwork_i );

    failed = compare_shseqr( h, h_i, wr, wr_i, wi, wi_i, z, z_i, info, info_i,
                             compz, ldh, ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to shseqr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to shseqr\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_shseqr( LAPACK_COL_MAJOR, job_i, compz_i, n_i, ilo_i,
                             ihi_i, h_i, ldh_i, wr_i, wi_i, z_i, ldz_i );

    failed = compare_shseqr( h, h_i, wr, wr_i, wi, wi_i, z, z_i, info, info_i,
                             compz, ldh, ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to shseqr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to shseqr\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, z_i, ldz, z_r, n+2 );
    }
    info_i = LAPACKE_shseqr_work( LAPACK_ROW_MAJOR, job_i, compz_i, n_i, ilo_i,
                                  ihi_i, h_r, ldh_r, wr_i, wi_i, z_r, ldz_r,
                                  work_i, lwork_i );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, n, h_r, n+2, h_i, ldh );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, n, z_r, n+2, z_i, ldz );
    }

    failed = compare_shseqr( h, h_i, wr, wr_i, wi, wi_i, z, z_i, info, info_i,
                             compz, ldh, ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to shseqr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to shseqr\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, n, z_i, ldz, z_r, n+2 );
    }
    info_i = LAPACKE_shseqr( LAPACK_ROW_MAJOR, job_i, compz_i, n_i, ilo_i,
                             ihi_i, h_r, ldh_r, wr_i, wi_i, z_r, ldz_r );

    LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, n, h_r, n+2, h_i, ldh );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_sge_trans( LAPACK_ROW_MAJOR, n, n, z_r, n+2, z_i, ldz );
    }

    failed = compare_shseqr( h, h_i, wr, wr_i, wi, wi_i, z, z_i, info, info_i,
                             compz, ldh, ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to shseqr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to shseqr\n" );
    }

    /* Release memory */
    if( h != NULL ) {
        LAPACKE_free( h );
    }
    if( h_i != NULL ) {
        LAPACKE_free( h_i );
    }
    if( h_r != NULL ) {
        LAPACKE_free( h_r );
    }
    if( h_save != NULL ) {
        LAPACKE_free( h_save );
    }
    if( wr != NULL ) {
        LAPACKE_free( wr );
    }
    if( wr_i != NULL ) {
        LAPACKE_free( wr_i );
    }
    if( wr_save != NULL ) {
        LAPACKE_free( wr_save );
    }
    if( wi != NULL ) {
        LAPACKE_free( wi );
    }
    if( wi_i != NULL ) {
        LAPACKE_free( wi_i );
    }
    if( wi_save != NULL ) {
        LAPACKE_free( wi_save );
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

/* Auxiliary function: shseqr scalar parameters initialization */
static void init_scalars_shseqr( char *job, char *compz, lapack_int *n,
                                 lapack_int *ilo, lapack_int *ihi,
                                 lapack_int *ldh, lapack_int *ldz,
                                 lapack_int *lwork )
{
    *job = 'S';
    *compz = 'I';
    *n = 4;
    *ilo = 1;
    *ihi = 4;
    *ldh = 8;
    *ldz = 8;
    *lwork = 8;

    return;
}

/* Auxiliary functions: shseqr array parameters initialization */
static void init_h( lapack_int size, float *h ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        h[i] = 0;
    }
    h[0] = 3.499999940e-001;  /* h[0,0] */
    h[8] = -1.159999967e-001;  /* h[0,1] */
    h[16] = -3.885999918e-001;  /* h[0,2] */
    h[24] = -2.942000031e-001;  /* h[0,3] */
    h[1] = -5.139999986e-001;  /* h[1,0] */
    h[9] = 1.225000024e-001;  /* h[1,1] */
    h[17] = 1.004000008e-001;  /* h[1,2] */
    h[25] = 1.125999987e-001;  /* h[1,3] */
    h[2] = 0.000000000e+000;  /* h[2,0] */
    h[10] = 6.442999840e-001;  /* h[2,1] */
    h[18] = -1.357000023e-001;  /* h[2,2] */
    h[26] = -9.769999981e-002;  /* h[2,3] */
    h[3] = 0.000000000e+000;  /* h[3,0] */
    h[11] = 0.000000000e+000;  /* h[3,1] */
    h[19] = 4.262000024e-001;  /* h[3,2] */
    h[27] = 1.632000059e-001;  /* h[3,3] */
}
static void init_wr( lapack_int size, float *wr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        wr[i] = 0;
    }
}
static void init_wi( lapack_int size, float *wi ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        wi[i] = 0;
    }
}
static void init_z( lapack_int size, float *z ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        z[i] = 0;
    }
    z[0] = 0.000000000e+000;  /* z[0,0] */
    z[8] = 0.000000000e+000;  /* z[0,1] */
    z[16] = 0.000000000e+000;  /* z[0,2] */
    z[24] = 0.000000000e+000;  /* z[0,3] */
    z[1] = 0.000000000e+000;  /* z[1,0] */
    z[9] = 0.000000000e+000;  /* z[1,1] */
    z[17] = 0.000000000e+000;  /* z[1,2] */
    z[25] = 0.000000000e+000;  /* z[1,3] */
    z[2] = 0.000000000e+000;  /* z[2,0] */
    z[10] = 0.000000000e+000;  /* z[2,1] */
    z[18] = 0.000000000e+000;  /* z[2,2] */
    z[26] = 0.000000000e+000;  /* z[2,3] */
    z[3] = 0.000000000e+000;  /* z[3,0] */
    z[11] = 0.000000000e+000;  /* z[3,1] */
    z[19] = 0.000000000e+000;  /* z[3,2] */
    z[27] = 0.000000000e+000;  /* z[3,3] */
}
static void init_work( lapack_int size, float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}

/* Auxiliary function: C interface to shseqr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_shseqr( float *h, float *h_i, float *wr, float *wr_i,
                           float *wi, float *wi_i, float *z, float *z_i,
                           lapack_int info, lapack_int info_i, char compz,
                           lapack_int ldh, lapack_int ldz, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldh*n; i++ ) {
        failed += compare_floats(h[i],h_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_floats(wr[i],wr_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_floats(wi[i],wi_i[i]);
    }
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        for( i = 0; i < ldz*n; i++ ) {
            failed += compare_floats(z[i],z_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
