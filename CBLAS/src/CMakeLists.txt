set(ALLBLAS cblas_globals.c cblas_xerbla.c xerbla.c)

#
#
# CBLAS routines
#
# Level 1
#
#

#
# All source files for single real precision
#
set(SLEV1 cblas_srotg.c cblas_srotmg.c cblas_srot.c cblas_srotm.c
          cblas_sswap.c cblas_sscal.c  cblas_scopy.c cblas_saxpy.c
          cblas_sdot.c  cblas_sdsdot.c cblas_snrm2.c cblas_sasum.c)

#
# All source files for double real precision
#
set(DLEV1 cblas_drotg.c cblas_drotmg.c cblas_drot.c cblas_drotm.c
          cblas_dswap.c cblas_dscal.c  cblas_dcopy.c cblas_daxpy.c
          cblas_ddot.c  cblas_dsdot.c cblas_dnrm2.c cblas_dasum.c
          cblas_idamax.c)

#
# All source files for single complex precision
#
set(CLEV1 cblas_cswap.c cblas_cscal.c cblas_csscal.c cblas_ccopy.c
          cblas_caxpy.c cblas_cdotu_sub.c cblas_cdotc_sub.c
          cblas_icamax.c)

#
# All source files for double complex precision
#
set(ZLEV1 cblas_zswap.c cblas_zscal.c cblas_zdscal.c cblas_zcopy.c
          cblas_zaxpy.c cblas_zdotu_sub.c cblas_zdotc_sub.c cblas_dznrm2.c
          cblas_dzasum.c cblas_izamax.c)

#
# Common files for single / complex precision
#
set(SCLEV1 cblas_scasum.c cblas_scnrm2.c)

#
# All source files
#
set(ALEV1 ${SLEV1} ${DLEV1} ${CLEV1} ${ZLEV1} ${SCLEV1})

#
#
# CBLAS routines
#
# Level 2
#
#

#
# All source files for single real precision
#
set(SLEV2 cblas_sgemv.c cblas_sgbmv.c cblas_sger.c cblas_ssbmv.c  cblas_sspmv.c
          cblas_sspr.c  cblas_sspr2.c cblas_ssymv.c cblas_ssyr.c  cblas_ssyr2.c
          cblas_stbmv.c cblas_stbsv.c cblas_stpmv.c cblas_stpsv.c  cblas_strmv.c
          cblas_strsv.c)

#
# All source files for double real precision
#
set(DLEV2 cblas_dgemv.c cblas_dgbmv.c cblas_dger.c cblas_dsbmv.c  cblas_dspmv.c
          cblas_dspr.c  cblas_dspr2.c cblas_dsymv.c cblas_dsyr.c  cblas_dsyr2.c
          cblas_dtbmv.c cblas_dtbsv.c cblas_dtpmv.c cblas_dtpsv.c  cblas_dtrmv.c
          cblas_dtrsv.c)

#
# All source files for single complex precision
#
set(CLEV2 cblas_cgemv.c cblas_cgbmv.c cblas_chemv.c cblas_chbmv.c cblas_chpmv.c
          cblas_ctrmv.c cblas_ctbmv.c cblas_ctpmv.c cblas_ctrsv.c cblas_ctbsv.c
          cblas_ctpsv.c cblas_cgeru.c cblas_cgerc.c cblas_cher.c  cblas_cher2.c
          cblas_chpr.c  cblas_chpr2.c)

#
# All source files for double complex precision
#
set(ZLEV2 cblas_zgemv.c cblas_zgbmv.c cblas_zhemv.c cblas_zhbmv.c cblas_zhpmv.c
          cblas_ztrmv.c cblas_ztbmv.c cblas_ztpmv.c cblas_ztrsv.c cblas_ztbsv.c
          cblas_ztpsv.c cblas_zgeru.c cblas_zgerc.c cblas_zher.c  cblas_zher2.c
          cblas_zhpr.c  cblas_zhpr2.c)

#
# All source files
#
set(ALEV2 ${SLEV2} ${DLEV2} ${CLEV2} ${ZLEV2})

#
#
# CBLAS routines
#
# Level 3
#
#

#
# All source files for single real precision
#
set(SLEV3 cblas_sgemm.c cblas_ssymm.c cblas_ssyrk.c cblas_ssyr2k.c cblas_strmm.c
          cblas_strsm.c)

#
# All source files for double real precision
#
set(DLEV3 cblas_dgemm.c cblas_dsymm.c cblas_dsyrk.c cblas_dsyr2k.c cblas_dtrmm.c
          cblas_dtrsm.c)

#
# All source files for single complex precision
#
set(CLEV3 cblas_cgemm.c cblas_csymm.c cblas_chemm.c cblas_cherk.c
          cblas_cher2k.c cblas_ctrmm.c cblas_ctrsm.c  cblas_csyrk.c
          cblas_csyr2k.c)

#
# All source files for double complex precision
#
set(ZLEV3 cblas_zgemm.c cblas_zsymm.c cblas_zhemm.c cblas_zherk.c
          cblas_zher2k.c cblas_ztrmm.c cblas_ztrsm.c  cblas_zsyrk.c
          cblas_zsyr2k.c)

#
# All source files
#
set(ALEV3 ${SLEV3} ${DLEV3} ${CLEV3} ${ZLEV3})

set(ALLOBJ ${ALEV1} ${ALEV2} ${ALEV3} ${ALLBLAS})

#add_library(cblas ${ALLOBJ})
add_library(cblas ${SLEV1} ${DLEV1} ${SLEV2} ${DLEV2} ${SLEV3} ${DLEV3} ${ALLBLAS})
target_link_libraries(cblas blas)
set_target_properties(cblas PROPERTIES COMPILE_FLAGS ${CBLAS_FAIL_FLAGS})

