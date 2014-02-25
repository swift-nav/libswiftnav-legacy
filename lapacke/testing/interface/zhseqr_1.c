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
* zhseqr_1 is the test program for the C interface to LAPACK
* routine zhseqr
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

static void init_scalars_zhseqr( char *job, char *compz, lapack_int *n,
                                 lapack_int *ilo, lapack_int *ihi,
                                 lapack_int *ldh, lapack_int *ldz,
                                 lapack_int *lwork );
static void init_h( lapack_int size, lapack_complex_double *h );
static void init_w( lapack_int size, lapack_complex_double *w );
static void init_z( lapack_int size, lapack_complex_double *z );
static void init_work( lapack_int size, lapack_complex_double *work );
static int compare_zhseqr( lapack_complex_double *h, lapack_complex_double *h_i,
                           lapack_complex_double *w, lapack_complex_double *w_i,
                           lapack_complex_double *z, lapack_complex_double *z_i,
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
    lapack_complex_double *h = NULL, *h_i = NULL;
    lapack_complex_double *w = NULL, *w_i = NULL;
    lapack_complex_double *z = NULL, *z_i = NULL;
    lapack_complex_double *work = NULL, *work_i = NULL;
    lapack_complex_double *h_save = NULL;
    lapack_complex_double *w_save = NULL;
    lapack_complex_double *z_save = NULL;
    lapack_complex_double *h_r = NULL;
    lapack_complex_double *z_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_zhseqr( &job, &compz, &n, &ilo, &ihi, &ldh, &ldz, &lwork );
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
    h = (lapack_complex_double *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_double) );
    w = (lapack_complex_double *)
        LAPACKE_malloc( n * sizeof(lapack_complex_double) );
    z = (lapack_complex_double *)
        LAPACKE_malloc( ldz*n * sizeof(lapack_complex_double) );
    work = (lapack_complex_double *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_double) );

    /* Allocate memory for the C interface function arrays */
    h_i = (lapack_complex_double *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_double) );
    w_i = (lapack_complex_double *)
        LAPACKE_malloc( n * sizeof(lapack_complex_double) );
    z_i = (lapack_complex_double *)
        LAPACKE_malloc( ldz*n * sizeof(lapack_complex_double) );
    work_i = (lapack_complex_double *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_double) );

    /* Allocate memory for the backup arrays */
    h_save = (lapack_complex_double *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_double) );
    w_save = (lapack_complex_double *)
        LAPACKE_malloc( n * sizeof(lapack_complex_double) );
    z_save = (lapack_complex_double *)
        LAPACKE_malloc( ldz*n * sizeof(lapack_complex_double) );

    /* Allocate memory for the row-major arrays */
    h_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );
    z_r = (lapack_complex_double *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_double) );

    /* Initialize input arrays */
    init_h( ldh*n, h );
    init_w( n, w );
    init_z( ldz*n, z );
    init_work( lwork, work );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldh*n; i++ ) {
        h_save[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        w_save[i] = w[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_save[i] = z[i];
    }

    /* Call the LAPACK routine */
    zhseqr_( &job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, work, &lwork,
             &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h_save[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_zhseqr_work( LAPACK_COL_MAJOR, job_i, compz_i, n_i, ilo_i,
                                  ihi_i, h_i, ldh_i, w_i, z_i, ldz_i, work_i,
                                  lwork_i );

    failed = compare_zhseqr( h, h_i, w, w_i, z, z_i, info, info_i, compz, ldh,
                             ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to zhseqr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to zhseqr\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h_save[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }
    info_i = LAPACKE_zhseqr( LAPACK_COL_MAJOR, job_i, compz_i, n_i, ilo_i,
                             ihi_i, h_i, ldh_i, w_i, z_i, ldz_i );

    failed = compare_zhseqr( h, h_i, w, w_i, z, z_i, info, info_i, compz, ldh,
                             ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to zhseqr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to zhseqr\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h_save[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, z_i, ldz, z_r, n+2 );
    }
    info_i = LAPACKE_zhseqr_work( LAPACK_ROW_MAJOR, job_i, compz_i, n_i, ilo_i,
                                  ihi_i, h_r, ldh_r, w_i, z_r, ldz_r, work_i,
                                  lwork_i );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, n, h_r, n+2, h_i, ldh );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, n, z_r, n+2, z_i, ldz );
    }

    failed = compare_zhseqr( h, h_i, w, w_i, z, z_i, info, info_i, compz, ldh,
                             ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to zhseqr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to zhseqr\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h_save[i];
    }
    for( i = 0; i < n; i++ ) {
        w_i[i] = w_save[i];
    }
    for( i = 0; i < ldz*n; i++ ) {
        z_i[i] = z_save[i];
    }
    for( i = 0; i < lwork; i++ ) {
        work_i[i] = work[i];
    }

    /* Init row_major arrays */
    LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, z_i, ldz, z_r, n+2 );
    }
    info_i = LAPACKE_zhseqr( LAPACK_ROW_MAJOR, job_i, compz_i, n_i, ilo_i,
                             ihi_i, h_r, ldh_r, w_i, z_r, ldz_r );

    LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, n, h_r, n+2, h_i, ldh );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_zge_trans( LAPACK_ROW_MAJOR, n, n, z_r, n+2, z_i, ldz );
    }

    failed = compare_zhseqr( h, h_i, w, w_i, z, z_i, info, info_i, compz, ldh,
                             ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to zhseqr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to zhseqr\n" );
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
    if( w != NULL ) {
        LAPACKE_free( w );
    }
    if( w_i != NULL ) {
        LAPACKE_free( w_i );
    }
    if( w_save != NULL ) {
        LAPACKE_free( w_save );
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

/* Auxiliary function: zhseqr scalar parameters initialization */
static void init_scalars_zhseqr( char *job, char *compz, lapack_int *n,
                                 lapack_int *ilo, lapack_int *ihi,
                                 lapack_int *ldh, lapack_int *ldz,
                                 lapack_int *lwork )
{
    *job = 'S';
    *compz = 'V';
    *n = 4;
    *ilo = 2;
    *ihi = 3;
    *ldh = 8;
    *ldz = 8;
    *lwork = 512;

    return;
}

/* Auxiliary functions: zhseqr array parameters initialization */
static void init_h( lapack_int size, lapack_complex_double *h ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        h[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    h[0] = lapack_make_complex_double( -1.25000000000000000e+000,
                                       7.50000000000000000e-001 );
    h[8] = lapack_make_complex_double( 1.38999999999999990e+000,
                                       3.97000000000000020e+000 );
    h[16] = lapack_make_complex_double( 8.28365051039430260e-001,
                                        7.39561433163036420e+000 );
    h[24] = lapack_make_complex_double( -2.08999999999999990e+000,
                                        7.55999999999999960e+000 );
    h[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    h[9] = lapack_make_complex_double( -2.50000000000000000e+000,
                                       -5.00000000000000000e-001 );
    h[17] = lapack_make_complex_double( 9.01376551729519800e-001,
                                        4.50688275864752350e-003 );
    h[25] = lapack_make_complex_double( -8.06000000000000050e+000,
                                        -1.24000000000000000e+000 );
    h[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    h[10] = lapack_make_complex_double( 1.10941425986869310e+000,
                                        0.00000000000000000e+000 );
    h[18] = lapack_make_complex_double( -2.50000000000000000e+000,
                                        -5.00000000000000000e-001 );
    h[26] = lapack_make_complex_double( -1.05960419162011950e+001,
                                        -4.66480393051061170e+000 );
    h[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    h[11] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    h[19] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    h[27] = lapack_make_complex_double( 1.50000000000000000e+000,
                                        -2.75000000000000000e+000 );
}
static void init_w( lapack_int size, lapack_complex_double *w ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        w[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
}
static void init_z( lapack_int size, lapack_complex_double *z ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        z[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
    z[0] = lapack_make_complex_double( 1.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    z[8] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    z[16] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    z[24] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    z[1] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    z[9] = lapack_make_complex_double( 1.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    z[17] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    z[25] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    z[2] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    z[10] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    z[18] = lapack_make_complex_double( -8.29266427591158320e-001,
                                        -5.58853462072302240e-001 );
    z[26] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    z[3] = lapack_make_complex_double( 0.00000000000000000e+000,
                                       0.00000000000000000e+000 );
    z[11] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    z[19] = lapack_make_complex_double( 0.00000000000000000e+000,
                                        0.00000000000000000e+000 );
    z[27] = lapack_make_complex_double( 1.00000000000000000e+000,
                                        0.00000000000000000e+000 );
}
static void init_work( lapack_int size, lapack_complex_double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_double( 0.0, 0.0 );
    }
}

/* Auxiliary function: C interface to zhseqr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_zhseqr( lapack_complex_double *h, lapack_complex_double *h_i,
                           lapack_complex_double *w, lapack_complex_double *w_i,
                           lapack_complex_double *z, lapack_complex_double *z_i,
                           lapack_int info, lapack_int info_i, char compz,
                           lapack_int ldh, lapack_int ldz, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldh*n; i++ ) {
        failed += compare_complex_doubles(h[i],h_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_complex_doubles(w[i],w_i[i]);
    }
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        for( i = 0; i < ldz*n; i++ ) {
            failed += compare_complex_doubles(z[i],z_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
