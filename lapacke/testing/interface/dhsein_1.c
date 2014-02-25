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
* dhsein_1 is the test program for the C interface to LAPACK
* routine dhsein
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

static void init_scalars_dhsein( char *job, char *eigsrc, char *initv,
                                 lapack_int *n, lapack_int *ldh,
                                 lapack_int *ldvl, lapack_int *ldvr,
                                 lapack_int *mm );
static void init_select( lapack_int size, lapack_int *select );
static void init_h( lapack_int size, double *h );
static void init_wr( lapack_int size, double *wr );
static void init_wi( lapack_int size, double *wi );
static void init_vl( lapack_int size, double *vl );
static void init_vr( lapack_int size, double *vr );
static void init_work( lapack_int size, double *work );
static void init_ifaill( lapack_int size, lapack_int *ifaill );
static void init_ifailr( lapack_int size, lapack_int *ifailr );
static int compare_dhsein( lapack_int *select, lapack_int *select_i, double *wr,
                           double *wr_i, double *vl, double *vl_i, double *vr,
                           double *vr_i, lapack_int m, lapack_int m_i,
                           lapack_int *ifaill, lapack_int *ifaill_i,
                           lapack_int *ifailr, lapack_int *ifailr_i,
                           lapack_int info, lapack_int info_i, char job,
                           lapack_int ldvl, lapack_int ldvr, lapack_int mm,
                           lapack_int n );

int main(void)
{
    /* Local scalars */
    char job, job_i;
    char eigsrc, eigsrc_i;
    char initv, initv_i;
    lapack_int n, n_i;
    lapack_int ldh, ldh_i;
    lapack_int ldh_r;
    lapack_int ldvl, ldvl_i;
    lapack_int ldvl_r;
    lapack_int ldvr, ldvr_i;
    lapack_int ldvr_r;
    lapack_int mm, mm_i;
    lapack_int m, m_i;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    lapack_int *select = NULL, *select_i = NULL;
    double *h = NULL, *h_i = NULL;
    double *wr = NULL, *wr_i = NULL;
    double *wi = NULL, *wi_i = NULL;
    double *vl = NULL, *vl_i = NULL;
    double *vr = NULL, *vr_i = NULL;
    double *work = NULL, *work_i = NULL;
    lapack_int *ifaill = NULL, *ifaill_i = NULL;
    lapack_int *ifailr = NULL, *ifailr_i = NULL;
    lapack_int *select_save = NULL;
    double *wr_save = NULL;
    double *vl_save = NULL;
    double *vr_save = NULL;
    lapack_int *ifaill_save = NULL;
    lapack_int *ifailr_save = NULL;
    double *h_r = NULL;
    double *vl_r = NULL;
    double *vr_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dhsein( &job, &eigsrc, &initv, &n, &ldh, &ldvl, &ldvr, &mm );
    ldh_r = n+2;
    ldvl_r = mm+2;
    ldvr_r = mm+2;
    job_i = job;
    eigsrc_i = eigsrc;
    initv_i = initv;
    n_i = n;
    ldh_i = ldh;
    ldvl_i = ldvl;
    ldvr_i = ldvr;
    mm_i = mm;

    /* Allocate memory for the LAPACK routine arrays */
    select = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    h = (double *)LAPACKE_malloc( ldh*n * sizeof(double) );
    wr = (double *)LAPACKE_malloc( n * sizeof(double) );
    wi = (double *)LAPACKE_malloc( n * sizeof(double) );
    vl = (double *)LAPACKE_malloc( ldvl*mm * sizeof(double) );
    vr = (double *)LAPACKE_malloc( ldvr*mm * sizeof(double) );
    work = (double *)LAPACKE_malloc( (((n+2)*n)) * sizeof(double) );
    ifaill = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the C interface function arrays */
    select_i = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    h_i = (double *)LAPACKE_malloc( ldh*n * sizeof(double) );
    wr_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    wi_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    vl_i = (double *)LAPACKE_malloc( ldvl*mm * sizeof(double) );
    vr_i = (double *)LAPACKE_malloc( ldvr*mm * sizeof(double) );
    work_i = (double *)LAPACKE_malloc( (((n+2)*n)) * sizeof(double) );
    ifaill_i = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr_i = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the backup arrays */
    select_save = (lapack_int *)LAPACKE_malloc( n * sizeof(lapack_int) );
    wr_save = (double *)LAPACKE_malloc( n * sizeof(double) );
    vl_save = (double *)LAPACKE_malloc( ldvl*mm * sizeof(double) );
    vr_save = (double *)LAPACKE_malloc( ldvr*mm * sizeof(double) );
    ifaill_save = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );
    ifailr_save = (lapack_int *)LAPACKE_malloc( mm * sizeof(lapack_int) );

    /* Allocate memory for the row-major arrays */
    h_r = (double *)LAPACKE_malloc( n*(n+2) * sizeof(double) );
    vl_r = (double *)LAPACKE_malloc( n*(mm+2) * sizeof(double) );
    vr_r = (double *)LAPACKE_malloc( n*(mm+2) * sizeof(double) );

    /* Initialize input arrays */
    init_select( n, select );
    init_h( ldh*n, h );
    init_wr( n, wr );
    init_wi( n, wi );
    init_vl( ldvl*mm, vl );
    init_vr( ldvr*mm, vr );
    init_work( ((n+2)*n), work );
    init_ifaill( mm, ifaill );
    init_ifailr( mm, ifailr );

    /* Backup the ouptut arrays */
    for( i = 0; i < n; i++ ) {
        select_save[i] = select[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_save[i] = wr[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_save[i] = vl[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_save[i] = vr[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_save[i] = ifaill[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_save[i] = ifailr[i];
    }

    /* Call the LAPACK routine */
    dhsein_( &job, &eigsrc, &initv, select, &n, h, &ldh, wr, wi, vl, &ldvl, vr,
             &ldvr, &mm, &m, work, ifaill, ifailr, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select_save[i];
    }
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl_save[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr_save[i];
    }
    for( i = 0; i < ((n+2)*n); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_i[i] = ifaill_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_i[i] = ifailr_save[i];
    }
    info_i = LAPACKE_dhsein_work( LAPACK_COL_MAJOR, job_i, eigsrc_i, initv_i,
                                  select_i, n_i, h_i, ldh_i, wr_i, wi_i, vl_i,
                                  ldvl_i, vr_i, ldvr_i, mm_i, &m_i, work_i,
                                  ifaill_i, ifailr_i );

    failed = compare_dhsein( select, select_i, wr, wr_i, vl, vl_i, vr, vr_i, m,
                             m_i, ifaill, ifaill_i, ifailr, ifailr_i, info,
                             info_i, job, ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dhsein\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dhsein\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select_save[i];
    }
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl_save[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr_save[i];
    }
    for( i = 0; i < ((n+2)*n); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_i[i] = ifaill_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_i[i] = ifailr_save[i];
    }
    info_i = LAPACKE_dhsein( LAPACK_COL_MAJOR, job_i, eigsrc_i, initv_i,
                             select_i, n_i, h_i, ldh_i, wr_i, wi_i, vl_i,
                             ldvl_i, vr_i, ldvr_i, mm_i, &m_i, ifaill_i,
                             ifailr_i );

    failed = compare_dhsein( select, select_i, wr, wr_i, vl, vl_i, vr, vr_i, m,
                             m_i, ifaill, ifaill_i, ifailr, ifailr_i, info,
                             info_i, job, ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dhsein\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dhsein\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select_save[i];
    }
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl_save[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr_save[i];
    }
    for( i = 0; i < ((n+2)*n); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_i[i] = ifaill_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_i[i] = ifailr_save[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_dhsein_work( LAPACK_ROW_MAJOR, job_i, eigsrc_i, initv_i,
                                  select_i, n_i, h_r, ldh_r, wr_i, wi_i, vl_r,
                                  ldvl_r, vr_r, ldvr_r, mm_i, &m_i, work_i,
                                  ifaill_i, ifailr_i );

    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, mm, vl_r, mm+2, vl_i, ldvl );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, mm, vr_r, mm+2, vr_i, ldvr );
    }

    failed = compare_dhsein( select, select_i, wr, wr_i, vl, vl_i, vr, vr_i, m,
                             m_i, ifaill, ifaill_i, ifailr, ifailr_i, info,
                             info_i, job, ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dhsein\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dhsein\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        select_i[i] = select_save[i];
    }
    for( i = 0; i < ldh*n; i++ ) {
        h_i[i] = h[i];
    }
    for( i = 0; i < n; i++ ) {
        wr_i[i] = wr_save[i];
    }
    for( i = 0; i < n; i++ ) {
        wi_i[i] = wi[i];
    }
    for( i = 0; i < ldvl*mm; i++ ) {
        vl_i[i] = vl_save[i];
    }
    for( i = 0; i < ldvr*mm; i++ ) {
        vr_i[i] = vr_save[i];
    }
    for( i = 0; i < ((n+2)*n); i++ ) {
        work_i[i] = work[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifaill_i[i] = ifaill_save[i];
    }
    for( i = 0; i < mm; i++ ) {
        ifailr_i[i] = ifailr_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, n, h_i, ldh, h_r, n+2 );
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, mm, vl_i, ldvl, vl_r, mm+2 );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, mm, vr_i, ldvr, vr_r, mm+2 );
    }
    info_i = LAPACKE_dhsein( LAPACK_ROW_MAJOR, job_i, eigsrc_i, initv_i,
                             select_i, n_i, h_r, ldh_r, wr_i, wi_i, vl_r,
                             ldvl_r, vr_r, ldvr_r, mm_i, &m_i, ifaill_i,
                             ifailr_i );

    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, mm, vl_r, mm+2, vl_i, ldvl );
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, mm, vr_r, mm+2, vr_i, ldvr );
    }

    failed = compare_dhsein( select, select_i, wr, wr_i, vl, vl_i, vr, vr_i, m,
                             m_i, ifaill, ifaill_i, ifailr, ifailr_i, info,
                             info_i, job, ldvl, ldvr, mm, n );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dhsein\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dhsein\n" );
    }

    /* Release memory */
    if( select != NULL ) {
        LAPACKE_free( select );
    }
    if( select_i != NULL ) {
        LAPACKE_free( select_i );
    }
    if( select_save != NULL ) {
        LAPACKE_free( select_save );
    }
    if( h != NULL ) {
        LAPACKE_free( h );
    }
    if( h_i != NULL ) {
        LAPACKE_free( h_i );
    }
    if( h_r != NULL ) {
        LAPACKE_free( h_r );
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
    if( vl != NULL ) {
        LAPACKE_free( vl );
    }
    if( vl_i != NULL ) {
        LAPACKE_free( vl_i );
    }
    if( vl_r != NULL ) {
        LAPACKE_free( vl_r );
    }
    if( vl_save != NULL ) {
        LAPACKE_free( vl_save );
    }
    if( vr != NULL ) {
        LAPACKE_free( vr );
    }
    if( vr_i != NULL ) {
        LAPACKE_free( vr_i );
    }
    if( vr_r != NULL ) {
        LAPACKE_free( vr_r );
    }
    if( vr_save != NULL ) {
        LAPACKE_free( vr_save );
    }
    if( work != NULL ) {
        LAPACKE_free( work );
    }
    if( work_i != NULL ) {
        LAPACKE_free( work_i );
    }
    if( ifaill != NULL ) {
        LAPACKE_free( ifaill );
    }
    if( ifaill_i != NULL ) {
        LAPACKE_free( ifaill_i );
    }
    if( ifaill_save != NULL ) {
        LAPACKE_free( ifaill_save );
    }
    if( ifailr != NULL ) {
        LAPACKE_free( ifailr );
    }
    if( ifailr_i != NULL ) {
        LAPACKE_free( ifailr_i );
    }
    if( ifailr_save != NULL ) {
        LAPACKE_free( ifailr_save );
    }

    return 0;
}

/* Auxiliary function: dhsein scalar parameters initialization */
static void init_scalars_dhsein( char *job, char *eigsrc, char *initv,
                                 lapack_int *n, lapack_int *ldh,
                                 lapack_int *ldvl, lapack_int *ldvr,
                                 lapack_int *mm )
{
    *job = 'R';
    *eigsrc = 'Q';
    *initv = 'N';
    *n = 4;
    *ldh = 8;
    *ldvl = 8;
    *ldvr = 8;
    *mm = 4;

    return;
}

/* Auxiliary functions: dhsein array parameters initialization */
static void init_select( lapack_int size, lapack_int *select ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        select[i] = 0;
    }
    select[0] = 0;
    select[1] = -1;
    select[2] = -1;
    select[3] = -1;
}
static void init_h( lapack_int size, double *h ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        h[i] = 0;
    }
    h[0] = 3.49999999999999980e-001;  /* h[0,0] */
    h[8] = -1.15952429620503390e-001;  /* h[0,1] */
    h[16] = -3.88601034323321160e-001;  /* h[0,2] */
    h[24] = -2.94184075347302120e-001;  /* h[0,3] */
    h[1] = -5.14003891035855980e-001;  /* h[1,0] */
    h[9] = 1.22486752460257420e-001;  /* h[1,1] */
    h[17] = 1.00359789682150170e-001;  /* h[1,2] */
    h[25] = 1.12561879970531830e-001;  /* h[1,3] */
    h[2] = -7.28472128292762870e-001;  /* h[2,0] */
    h[10] = 6.44263618527061930e-001;  /* h[2,1] */
    h[18] = -1.35700171757113630e-001;  /* h[2,2] */
    h[26] = -9.76816227049334410e-002;  /* h[2,3] */
    h[3] = 4.13904618348160720e-001;  /* h[3,0] */
    h[11] = -1.66544579490569860e-001;  /* h[3,1] */
    h[19] = 4.26244372207844720e-001;  /* h[3,2] */
    h[27] = 1.63213419296856090e-001;  /* h[3,3] */
}
static void init_wr( lapack_int size, double *wr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        wr[i] = 0;
    }
    wr[0] = 7.99482122586209880e-001;
    wr[1] = -9.94124532950746570e-002;
    wr[2] = -9.94124532950746570e-002;
    wr[3] = -1.00657215996058770e-001;
}
static void init_wi( lapack_int size, double *wi ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        wi[i] = 0;
    }
    wi[0] = 0.00000000000000000e+000;
    wi[1] = 4.00792471989754480e-001;
    wi[2] = -4.00792471989754480e-001;
    wi[3] = 0.00000000000000000e+000;
}
static void init_vl( lapack_int size, double *vl ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vl[i] = 0;
    }
    vl[0] = 0.00000000000000000e+000;  /* vl[0,0] */
    vl[8] = 0.00000000000000000e+000;  /* vl[0,1] */
    vl[16] = 0.00000000000000000e+000;  /* vl[0,2] */
    vl[24] = 0.00000000000000000e+000;  /* vl[0,3] */
    vl[1] = 0.00000000000000000e+000;  /* vl[1,0] */
    vl[9] = 0.00000000000000000e+000;  /* vl[1,1] */
    vl[17] = 0.00000000000000000e+000;  /* vl[1,2] */
    vl[25] = 0.00000000000000000e+000;  /* vl[1,3] */
    vl[2] = 0.00000000000000000e+000;  /* vl[2,0] */
    vl[10] = 0.00000000000000000e+000;  /* vl[2,1] */
    vl[18] = 0.00000000000000000e+000;  /* vl[2,2] */
    vl[26] = 0.00000000000000000e+000;  /* vl[2,3] */
    vl[3] = 0.00000000000000000e+000;  /* vl[3,0] */
    vl[11] = 0.00000000000000000e+000;  /* vl[3,1] */
    vl[19] = 0.00000000000000000e+000;  /* vl[3,2] */
    vl[27] = 0.00000000000000000e+000;  /* vl[3,3] */
}
static void init_vr( lapack_int size, double *vr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        vr[i] = 0;
    }
    vr[0] = 0.00000000000000000e+000;  /* vr[0,0] */
    vr[8] = 0.00000000000000000e+000;  /* vr[0,1] */
    vr[16] = 0.00000000000000000e+000;  /* vr[0,2] */
    vr[24] = 0.00000000000000000e+000;  /* vr[0,3] */
    vr[1] = 0.00000000000000000e+000;  /* vr[1,0] */
    vr[9] = 0.00000000000000000e+000;  /* vr[1,1] */
    vr[17] = 0.00000000000000000e+000;  /* vr[1,2] */
    vr[25] = 0.00000000000000000e+000;  /* vr[1,3] */
    vr[2] = 0.00000000000000000e+000;  /* vr[2,0] */
    vr[10] = 0.00000000000000000e+000;  /* vr[2,1] */
    vr[18] = 0.00000000000000000e+000;  /* vr[2,2] */
    vr[26] = 0.00000000000000000e+000;  /* vr[2,3] */
    vr[3] = 0.00000000000000000e+000;  /* vr[3,0] */
    vr[11] = 0.00000000000000000e+000;  /* vr[3,1] */
    vr[19] = 0.00000000000000000e+000;  /* vr[3,2] */
    vr[27] = 0.00000000000000000e+000;  /* vr[3,3] */
}
static void init_work( lapack_int size, double *work ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        work[i] = 0;
    }
}
static void init_ifaill( lapack_int size, lapack_int *ifaill ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ifaill[i] = 0;
    }
}
static void init_ifailr( lapack_int size, lapack_int *ifailr ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        ifailr[i] = 0;
    }
}

/* Auxiliary function: C interface to dhsein results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dhsein( lapack_int *select, lapack_int *select_i, double *wr,
                           double *wr_i, double *vl, double *vl_i, double *vr,
                           double *vr_i, lapack_int m, lapack_int m_i,
                           lapack_int *ifaill, lapack_int *ifaill_i,
                           lapack_int *ifailr, lapack_int *ifailr_i,
                           lapack_int info, lapack_int info_i, char job,
                           lapack_int ldvl, lapack_int ldvr, lapack_int mm,
                           lapack_int n )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < n; i++ ) {
        failed += (select[i] == select_i[i]) ? 0 : 1;
    }
    for( i = 0; i < n; i++ ) {
        failed += compare_doubles(wr[i],wr_i[i]);
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        for( i = 0; i < ldvl*mm; i++ ) {
            failed += compare_doubles(vl[i],vl_i[i]);
        }
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        for( i = 0; i < ldvr*mm; i++ ) {
            failed += compare_doubles(vr[i],vr_i[i]);
        }
    }
    failed += (m == m_i) ? 0 : 1;
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'l' ) ) {
        for( i = 0; i < mm; i++ ) {
            failed += (ifaill[i] == ifaill_i[i]) ? 0 : 1;
        }
    }
    if( LAPACKE_lsame( job, 'b' ) || LAPACKE_lsame( job, 'r' ) ) {
        for( i = 0; i < mm; i++ ) {
            failed += (ifailr[i] == ifailr_i[i]) ? 0 : 1;
        }
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
