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
* chseqr_3 is the test program for the C interface to LAPACK
* routine chseqr
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

static void init_scalars_chseqr( char *job, char *compz, lapack_int *n,
                                 lapack_int *ilo, lapack_int *ihi,
                                 lapack_int *ldh, lapack_int *ldz,
                                 lapack_int *lwork );
static void init_h( lapack_int size, lapack_complex_float *h );
static void init_w( lapack_int size, lapack_complex_float *w );
static void init_z( lapack_int size, lapack_complex_float *z );
static void init_work( lapack_int size, lapack_complex_float *work );
static int compare_chseqr( lapack_complex_float *h, lapack_complex_float *h_i,
                           lapack_complex_float *w, lapack_complex_float *w_i,
                           lapack_complex_float *z, lapack_complex_float *z_i,
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
    lapack_complex_float *h = NULL, *h_i = NULL;
    lapack_complex_float *w = NULL, *w_i = NULL;
    lapack_complex_float *z = NULL, *z_i = NULL;
    lapack_complex_float *work = NULL, *work_i = NULL;
    lapack_complex_float *h_save = NULL;
    lapack_complex_float *w_save = NULL;
    lapack_complex_float *z_save = NULL;
    lapack_complex_float *h_r = NULL;
    lapack_complex_float *z_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_chseqr( &job, &compz, &n, &ilo, &ihi, &ldh, &ldz, &lwork );
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
    h = (lapack_complex_float *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_float) );
    w = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );
    z = (lapack_complex_float *)
        LAPACKE_malloc( ldz*n * sizeof(lapack_complex_float) );
    work = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the C interface function arrays */
    h_i = (lapack_complex_float *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_float) );
    w_i = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );
    z_i = (lapack_complex_float *)
        LAPACKE_malloc( ldz*n * sizeof(lapack_complex_float) );
    work_i = (lapack_complex_float *)
        LAPACKE_malloc( lwork * sizeof(lapack_complex_float) );

    /* Allocate memory for the backup arrays */
    h_save = (lapack_complex_float *)
        LAPACKE_malloc( ldh*n * sizeof(lapack_complex_float) );
    w_save = (lapack_complex_float *)
        LAPACKE_malloc( n * sizeof(lapack_complex_float) );
    z_save = (lapack_complex_float *)
        LAPACKE_malloc( ldz*n * sizeof(lapack_complex_float) );

    /* Allocate memory for the row-major arrays */
    h_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );
    z_r = (lapack_complex_float *)
        LAPACKE_malloc( n*(n+2) * sizeof(lapack_complex_float) );

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
    chseqr_( &job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, work, &lwork,
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
    info_i = LAPACKE_chseqr_work( LAPACK_COL_MAJOR, job_i, compz_i, n_i, ilo_i,
                                  ihi_i, h_i, ldh_i, w_i, z_i, ldz_i, work_i,
                                  lwork_i );

    failed = compare_chseqr( h, h_i, w, w_i, z, z_i, info, info_i, compz, ldh,
                             ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to chseqr\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to chseqr\n" );
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
    info_i = LAPACKE_chseqr( LAPACK_COL_MAJOR, job_i, compz_i, n_i, ilo_i,
                             ihi_i, h_i, ldh_i, w_i, z_i, ldz_i );

    failed = compare_chseqr( h, h_i, w, w_i, z, z_i, info, info_i, compz, ldh,
                             ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to chseqr\n" );
    } else {
        printf( "FAILED: column-major high-level interface to chseqr\n" );
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

    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, z_i, ldz, z_r, n+2 );
    }
    info_i = LAPACKE_chseqr_work( LAPACK_ROW_MAJOR, job_i, compz_i, n_i, ilo_i,
                                  ihi_i, h_r, ldh_r, w_i, z_r, ldz_r, work_i,
                                  lwork_i );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, h_r, n+2, h_i, ldh );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, z_r, n+2, z_i, ldz );
    }

    failed = compare_chseqr( h, h_i, w, w_i, z, z_i, info, info_i, compz, ldh,
                             ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to chseqr\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to chseqr\n" );
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
    LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, z_i, ldz, z_r, n+2 );
    }
    info_i = LAPACKE_chseqr( LAPACK_ROW_MAJOR, job_i, compz_i, n_i, ilo_i,
                             ihi_i, h_r, ldh_r, w_i, z_r, ldz_r );

    LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, h_r, n+2, h_i, ldh );
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        LAPACKE_cge_trans( LAPACK_ROW_MAJOR, n, n, z_r, n+2, z_i, ldz );
    }

    failed = compare_chseqr( h, h_i, w, w_i, z, z_i, info, info_i, compz, ldh,
                             ldz, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to chseqr\n" );
    } else {
        printf( "FAILED: row-major high-level interface to chseqr\n" );
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

/* Auxiliary function: chseqr scalar parameters initialization */
static void init_scalars_chseqr( char *job, char *compz, lapack_int *n,
                                 lapack_int *ilo, lapack_int *ihi,
                                 lapack_int *ldh, lapack_int *ldz,
                                 lapack_int *lwork )
{
    *job = 'S';
    *compz = 'V';
    *n = 4;
    *ilo = 1;
    *ihi = 4;
    *ldh = 8;
    *ldz = 8;
    *lwork = 448;

    return;
}

/* Auxiliary functions: chseqr array parameters initialization */
static void init_h( lapack_int size, lapack_complex_float *h ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        h[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    h[0] = lapack_make_complex_float( -3.970000029e+000, -5.039999962e+000 );
    h[8] = lapack_make_complex_float( -1.131804943e+000, -2.569304705e+000 );
    h[16] = lapack_make_complex_float( -4.602741241e+000, -1.426316500e-001 );
    h[24] = lapack_make_complex_float( -1.424912333e+000, 1.732983828e+000 );
    h[1] = lapack_make_complex_float( -5.479653358e+000, 0.000000000e+000 );
    h[9] = lapack_make_complex_float( 1.858472466e+000, -1.550180435e+000 );
    h[17] = lapack_make_complex_float( 4.414464474e+000, -7.638237476e-001 );
    h[25] = lapack_make_complex_float( -4.805260897e-001, -1.197600603e+000 );
    h[2] = lapack_make_complex_float( 6.932221651e-001, -4.828752279e-001 );
    h[10] = lapack_make_complex_float( 6.267275810e+000, 0.000000000e+000 );
    h[18] = lapack_make_complex_float( -4.503800869e-001, -2.898204327e-002 );
    h[26] = lapack_make_complex_float( -1.346683741e+000, 1.657924891e+000 );
    h[3] = lapack_make_complex_float( -2.112946808e-001, 8.644121885e-002 );
    h[11] = lapack_make_complex_float( 1.242147088e-001, -2.289276123e-001 );
    h[19] = lapack_make_complex_float( -3.499985933e+000, 0.000000000e+000 );
    h[27] = lapack_make_complex_float( 2.561908484e+000, -3.370837450e+000 );
}
static void init_w( lapack_int size, lapack_complex_float *w ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        w[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}
static void init_z( lapack_int size, lapack_complex_float *z ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        z[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
    z[0] = lapack_make_complex_float( 1.000000000e+000, 0.000000000e+000 );
    z[8] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    z[16] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    z[24] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    z[1] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    z[9] = lapack_make_complex_float( -6.204771996e-002, 2.737399340e-001 );
    z[17] = lapack_make_complex_float( 4.994977713e-001, 6.440132260e-001 );
    z[25] = lapack_make_complex_float( 4.237119257e-001, -2.782686949e-001 );
    z[2] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    z[10] = lapack_make_complex_float( -6.040527821e-001, 7.025991082e-001 );
    z[18] = lapack_make_complex_float( -1.486803293e-001, -1.426575333e-001 );
    z[26] = lapack_make_complex_float( -2.633657455e-001, -1.722092628e-001 );
    z[3] = lapack_make_complex_float( 0.000000000e+000, 0.000000000e+000 );
    z[11] = lapack_make_complex_float( 2.007426172e-001, -1.496444941e-001 );
    z[19] = lapack_make_complex_float( -4.651779532e-001, 2.773107588e-001 );
    z[27] = lapack_make_complex_float( -2.482300401e-001, -7.631506920e-001 );
}
static void init_work( lapack_int size, lapack_complex_float *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = lapack_make_complex_float( 0.0f, 0.0f );
    }
}

/* Auxiliary function: C interface to chseqr results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_chseqr( lapack_complex_float *h, lapack_complex_float *h_i,
                           lapack_complex_float *w, lapack_complex_float *w_i,
                           lapack_complex_float *z, lapack_complex_float *z_i,
                           lapack_int info, lapack_int info_i, char compz,
                           lapack_int ldh, lapack_int ldz, lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldh*n; i++ ) {
        failed += compare_complex_floats(h[i],h_i[i]);
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_complex_floats(w[i],w_i[i]);
    }
    if( LAPACKE_lsame( compz, 'i' ) || LAPACKE_lsame( compz, 'v' ) ) {
        for( i = 0; i < ldz*n; i++ ) {
            failed += compare_complex_floats(z[i],z_i[i]);
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
