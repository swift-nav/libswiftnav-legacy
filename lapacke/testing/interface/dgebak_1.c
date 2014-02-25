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
* dgebak_1 is the test program for the C interface to LAPACK
* routine dgebak
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

static void init_scalars_dgebak( char *job, char *side, lapack_int *n,
                                 lapack_int *ilo, lapack_int *ihi,
                                 lapack_int *m, lapack_int *ldv );
static void init_scale( lapack_int size, double *scale );
static void init_v( lapack_int size, double *v );
static int compare_dgebak( double *v, double *v_i, lapack_int info,
                           lapack_int info_i, lapack_int ldv, lapack_int m );

int main(void)
{
    /* Local scalars */
    char job, job_i;
    char side, side_i;
    lapack_int n, n_i;
    lapack_int ilo, ilo_i;
    lapack_int ihi, ihi_i;
    lapack_int m, m_i;
    lapack_int ldv, ldv_i;
    lapack_int ldv_r;
    lapack_int info, info_i;
    lapack_int i;
    int failed;

    /* Local arrays */
    double *scale = NULL, *scale_i = NULL;
    double *v = NULL, *v_i = NULL;
    double *v_save = NULL;
    double *v_r = NULL;

    /* Iniitialize the scalar parameters */
    init_scalars_dgebak( &job, &side, &n, &ilo, &ihi, &m, &ldv );
    ldv_r = m+2;
    job_i = job;
    side_i = side;
    n_i = n;
    ilo_i = ilo;
    ihi_i = ihi;
    m_i = m;
    ldv_i = ldv;

    /* Allocate memory for the LAPACK routine arrays */
    scale = (double *)LAPACKE_malloc( n * sizeof(double) );
    v = (double *)LAPACKE_malloc( ldv*m * sizeof(double) );

    /* Allocate memory for the C interface function arrays */
    scale_i = (double *)LAPACKE_malloc( n * sizeof(double) );
    v_i = (double *)LAPACKE_malloc( ldv*m * sizeof(double) );

    /* Allocate memory for the backup arrays */
    v_save = (double *)LAPACKE_malloc( ldv*m * sizeof(double) );

    /* Allocate memory for the row-major arrays */
    v_r = (double *)LAPACKE_malloc( n*(m+2) * sizeof(double) );

    /* Initialize input arrays */
    init_scale( n, scale );
    init_v( ldv*m, v );

    /* Backup the ouptut arrays */
    for( i = 0; i < ldv*m; i++ ) {
        v_save[i] = v[i];
    }

    /* Call the LAPACK routine */
    dgebak_( &job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info );

    /* Initialize input data, call the column-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        scale_i[i] = scale[i];
    }
    for( i = 0; i < ldv*m; i++ ) {
        v_i[i] = v_save[i];
    }
    info_i = LAPACKE_dgebak_work( LAPACK_COL_MAJOR, job_i, side_i, n_i, ilo_i,
                                  ihi_i, scale_i, m_i, v_i, ldv_i );

    failed = compare_dgebak( v, v_i, info, info_i, ldv, m );
    if( failed == 0 ) {
        printf( "PASSED: column-major middle-level interface to dgebak\n" );
    } else {
        printf( "FAILED: column-major middle-level interface to dgebak\n" );
    }

    /* Initialize input data, call the column-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        scale_i[i] = scale[i];
    }
    for( i = 0; i < ldv*m; i++ ) {
        v_i[i] = v_save[i];
    }
    info_i = LAPACKE_dgebak( LAPACK_COL_MAJOR, job_i, side_i, n_i, ilo_i, ihi_i,
                             scale_i, m_i, v_i, ldv_i );

    failed = compare_dgebak( v, v_i, info, info_i, ldv, m );
    if( failed == 0 ) {
        printf( "PASSED: column-major high-level interface to dgebak\n" );
    } else {
        printf( "FAILED: column-major high-level interface to dgebak\n" );
    }

    /* Initialize input data, call the row-major middle-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        scale_i[i] = scale[i];
    }
    for( i = 0; i < ldv*m; i++ ) {
        v_i[i] = v_save[i];
    }

    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, m, v_i, ldv, v_r, m+2 );
    info_i = LAPACKE_dgebak_work( LAPACK_ROW_MAJOR, job_i, side_i, n_i, ilo_i,
                                  ihi_i, scale_i, m_i, v_r, ldv_r );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, m, v_r, m+2, v_i, ldv );

    failed = compare_dgebak( v, v_i, info, info_i, ldv, m );
    if( failed == 0 ) {
        printf( "PASSED: row-major middle-level interface to dgebak\n" );
    } else {
        printf( "FAILED: row-major middle-level interface to dgebak\n" );
    }

    /* Initialize input data, call the row-major high-level
     * interface to LAPACK routine and check the results */
    for( i = 0; i < n; i++ ) {
        scale_i[i] = scale[i];
    }
    for( i = 0; i < ldv*m; i++ ) {
        v_i[i] = v_save[i];
    }

    /* Init row_major arrays */
    LAPACKE_dge_trans( LAPACK_COL_MAJOR, n, m, v_i, ldv, v_r, m+2 );
    info_i = LAPACKE_dgebak( LAPACK_ROW_MAJOR, job_i, side_i, n_i, ilo_i, ihi_i,
                             scale_i, m_i, v_r, ldv_r );

    LAPACKE_dge_trans( LAPACK_ROW_MAJOR, n, m, v_r, m+2, v_i, ldv );

    failed = compare_dgebak( v, v_i, info, info_i, ldv, m );
    if( failed == 0 ) {
        printf( "PASSED: row-major high-level interface to dgebak\n" );
    } else {
        printf( "FAILED: row-major high-level interface to dgebak\n" );
    }

    /* Release memory */
    if( scale != NULL ) {
        LAPACKE_free( scale );
    }
    if( scale_i != NULL ) {
        LAPACKE_free( scale_i );
    }
    if( v != NULL ) {
        LAPACKE_free( v );
    }
    if( v_i != NULL ) {
        LAPACKE_free( v_i );
    }
    if( v_r != NULL ) {
        LAPACKE_free( v_r );
    }
    if( v_save != NULL ) {
        LAPACKE_free( v_save );
    }

    return 0;
}

/* Auxiliary function: dgebak scalar parameters initialization */
static void init_scalars_dgebak( char *job, char *side, lapack_int *n,
                                 lapack_int *ilo, lapack_int *ihi,
                                 lapack_int *m, lapack_int *ldv )
{
    *job = 'B';
    *side = 'R';
    *n = 4;
    *ilo = 2;
    *ihi = 4;
    *m = 4;
    *ldv = 8;

    return;
}

/* Auxiliary functions: dgebak array parameters initialization */
static void init_scale( lapack_int size, double *scale ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        scale[i] = 0;
    }
    scale[0] = 3.00000000000000000e+000;
    scale[1] = 4.00000000000000000e+000;
    scale[2] = 4.00000000000000000e+000;
    scale[3] = 5.00000000000000000e-001;
}
static void init_v( lapack_int size, double *v ) {
    lapack_int i;
    for( i = 0; i < size; i++ ) {
        v[i] = 0;
    }
    v[0] = 1.00000000000000000e+000;  /* v[0,0] */
    v[8] = -2.15709029223295460e-001;  /* v[0,1] */
    v[16] = -1.00000000000000000e+000;  /* v[0,2] */
    v[24] = -1.00000000000000000e+000;  /* v[0,3] */
    v[1] = 0.00000000000000000e+000;  /* v[1,0] */
    v[9] = 1.00000000000000000e+000;  /* v[1,1] */
    v[17] = -4.95308301949082260e-001;  /* v[1,2] */
    v[25] = 1.71833107632870240e-001;  /* v[1,3] */
    v[2] = 0.00000000000000000e+000;  /* v[2,0] */
    v[10] = -4.91006373625574700e-001;  /* v[2,1] */
    v[18] = -2.92210843245896300e-001;  /* v[2,2] */
    v[26] = -9.53731439243778790e-001;  /* v[2,3] */
    v[3] = 0.00000000000000000e+000;  /* v[3,0] */
    v[11] = -8.75122446424736690e-001;  /* v[3,1] */
    v[19] = -2.61487853687739190e-001;  /* v[3,2] */
    v[27] = 4.72487441336102100e-001;  /* v[3,3] */
}

/* Auxiliary function: C interface to dgebak results check */
/* Return value: 0 - test is passed, non-zero - test is failed */
static int compare_dgebak( double *v, double *v_i, lapack_int info,
                           lapack_int info_i, lapack_int ldv, lapack_int m )
{
    lapack_int i;
    int failed = 0;
    for( i = 0; i < ldv*m; i++ ) {
        failed += compare_doubles(v[i],v_i[i]);
    }
    failed += (info == info_i) ? 0 : 1;
    if( info != 0 || info_i != 0 ) {
        printf( "info=%d, info_i=%d\n",(int)info,(int)info_i );
    }

    return failed;
}
