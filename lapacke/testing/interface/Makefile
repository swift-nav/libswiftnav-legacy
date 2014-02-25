##############################################################################
# Copyright (c) 2010, Intel Corp.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   * Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of Intel Corporation nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.
##############################################################################
# Contents: Native C interface to LAPACK
# Author: Intel Corporation
# September, 2010
##############################################################################
# makefile for LAPACKE, used to test lapacke binary.
#
include ../../make.inc

TST = cbdsqr_1 \
      cbdsqr_2 \
      cgbbrd_1 \
      cgbcon_1 \
      cgbrfs_1 \
      cgbtrf_1 \
      cgbtrs_1 \
      cgebak_1 \
      cgebal_1 \
      cgebrd_1 \
      cgebrd_3 \
      cgebrd_4 \
      cgebrd_5 \
      cgecon_1 \
      cgehrd_1 \
      cgehrd_2 \
      cgehrd_3 \
      cgelqf_1 \
      cgelqf_3 \
      cgeqpf_1 \
      cgeqrf_1 \
      cgeqrf_3 \
      cgerfs_1 \
      cgetrf_1 \
      cgetri_1 \
      cgetrs_1 \
      chbtrd_1 \
      checon_1 \
      chegst_1 \
      cherfs_1 \
      chetrd_1 \
      chetrd_2 \
      chetrd_3 \
      chetrf_1 \
      chetri_1 \
      chetrs_1 \
      chpcon_1 \
      chpgst_1 \
      chprfs_1 \
      chptrd_1 \
      chptrd_2 \
      chptrf_1 \
      chptri_1 \
      chptrs_1 \
      chsein_1 \
      chseqr_1 \
      chseqr_2 \
      chseqr_3 \
      chseqr_4 \
      cpbcon_1 \
      cpbrfs_1 \
      cpbtrf_1 \
      cpbtrs_1 \
      cpocon_1 \
      cporfs_1 \
      cpotrf_1 \
      cpotri_1 \
      cpotrs_1 \
      cppcon_1 \
      cpprfs_1 \
      cpptrf_1 \
      cpptri_1 \
      cpptrs_1 \
      cpteqr_1 \
      cspcon_1 \
      csprfs_1 \
      csptrf_1 \
      csptri_1 \
      csptrs_1 \
      cstein_1 \
      csteqr_1 \
      csteqr_2 \
      csycon_1 \
      csyrfs_1 \
      csytrf_1 \
      csytri_1 \
      csytrs_1 \
      ctbcon_1 \
      ctbrfs_1 \
      ctbtrs_1 \
      ctpcon_1 \
      ctprfs_1 \
      ctptri_1 \
      ctptrs_1 \
      ctrcon_1 \
      ctrevc_1 \
      ctrevc_2 \
      ctrexc_1 \
      ctrrfs_1 \
      ctrsen_1 \
      ctrsna_1 \
      ctrsyl_1 \
      ctrtri_1 \
      ctrtrs_1 \
      cungbr_1 \
      cungbr_2 \
      cungbr_3 \
      cungbr_4 \
      cunghr_1 \
      cunghr_2 \
      cunglq_1 \
      cunglq_2 \
      cungqr_1 \
      cungqr_2 \
      cungtr_1 \
      cungtr_2 \
      cunmbr_1 \
      cunmbr_2 \
      cunmhr_1 \
      cunmlq_1 \
      cunmqr_1 \
      cunmqr_2 \
      cunmtr_1 \
      cupgtr_1 \
      cupmtr_1 \
      dbdsqr_1 \
      dbdsqr_2 \
      dbdsqr_3 \
      dgbbrd_1 \
      dgbcon_1 \
      dgbrfs_1 \
      dgbtrf_1 \
      dgbtrs_1 \
      dgebak_1 \
      dgebal_1 \
      dgebrd_1 \
      dgebrd_3 \
      dgebrd_4 \
      dgebrd_5 \
      dgecon_1 \
      dgehrd_1 \
      dgehrd_2 \
      dgehrd_3 \
      dgelqf_1 \
      dgelqf_3 \
      dgeqpf_1 \
      dgeqrf_1 \
      dgeqrf_2 \
      dgerfs_1 \
      dgetrf_1 \
      dgetri_1 \
      dgetrs_1 \
      dhsein_1 \
      dhseqr_1 \
      dhseqr_2 \
      dhseqr_3 \
      dhseqr_4 \
      dopgtr_1 \
      dopmtr_1 \
      dorgbr_1 \
      dorgbr_2 \
      dorgbr_3 \
      dorgbr_4 \
      dorghr_1 \
      dorghr_2 \
      dorglq_1 \
      dorglq_2 \
      dorgqr_1 \
      dorgqr_2 \
      dorgtr_1 \
      dormbr_1 \
      dormbr_2 \
      dormhr_1 \
      dormlq_1 \
      dormqr_1 \
      dormqr_2 \
      dormtr_1 \
      dpbcon_1 \
      dpbrfs_1 \
      dpbtrf_1 \
      dpbtrs_1 \
      dpocon_1 \
      dporfs_1 \
      dpotrf_1 \
      dpotrf_6 \
      dpotri_1 \
      dpotrs_1 \
      dppcon_1 \
      dpprfs_1 \
      dpptrf_1 \
      dpptrf_6 \
      dpptri_1 \
      dpptrs_1 \
      dpteqr_1 \
      dsbtrd_1 \
      dspcon_1 \
      dspgst_1 \
      dsprfs_1 \
      dsptrd_1 \
      dsptrd_3 \
      dsptrf_1 \
      dsptri_1 \
      dsptrs_1 \
      dstebz_1 \
      dstebz_2 \
      dstebz_3 \
      dstebz_4 \
      dstein_1 \
      dsteqr_1 \
      dsteqr_3 \
      dsteqr_4 \
      dsterf_1 \
      dsterf_2 \
      dsterf_3 \
      dsterf_4 \
      dsterf_5 \
      dsycon_1 \
      dsygst_1 \
      dsyrfs_1 \
      dsytrd_1 \
      dsytrd_3 \
      dsytrf_1 \
      dsytri_1 \
      dsytrs_1 \
      dtbcon_1 \
      dtbrfs_1 \
      dtbtrs_1 \
      dtpcon_1 \
      dtprfs_1 \
      dtptri_1 \
      dtptrs_1 \
      dtrcon_1 \
      dtrevc_1 \
      dtrevc_2 \
      dtrexc_1 \
      dtrrfs_1 \
      dtrsen_1 \
      dtrsna_1 \
      dtrsyl_1 \
      dtrtri_1 \
      dtrtrs_1 \
      sbdsqr_1 \
      sbdsqr_2 \
      sbdsqr_3 \
      sgbbrd_1 \
      sgbcon_1 \
      sgbrfs_1 \
      sgbtrf_1 \
      sgbtrs_1 \
      sgebak_1 \
      sgebal_1 \
      sgebrd_1 \
      sgebrd_3 \
      sgebrd_4 \
      sgebrd_5 \
      sgecon_1 \
      sgehrd_1 \
      sgehrd_2 \
      sgehrd_3 \
      sgelqf_1 \
      sgelqf_3 \
      sgeqpf_1 \
      sgeqrf_1 \
      sgeqrf_2 \
      sgerfs_1 \
      sgetrf_1 \
      sgetri_1 \
      sgetrs_1 \
      shsein_1 \
      shseqr_1 \
      shseqr_2 \
      shseqr_3 \
      shseqr_4 \
      sopgtr_1 \
      sopmtr_1 \
      sorgbr_1 \
      sorgbr_2 \
      sorgbr_3 \
      sorgbr_4 \
      sorghr_1 \
      sorghr_2 \
      sorglq_1 \
      sorglq_2 \
      sorgqr_1 \
      sorgqr_2 \
      sorgtr_1 \
      sormbr_1 \
      sormbr_2 \
      sormhr_1 \
      sormlq_1 \
      sormqr_1 \
      sormqr_2 \
      sormtr_1 \
      spbcon_1 \
      spbrfs_1 \
      spbtrf_1 \
      spbtrs_1 \
      spocon_1 \
      sporfs_1 \
      spotrf_1 \
      spotrf_6 \
      spotri_1 \
      spotrs_1 \
      sppcon_1 \
      spprfs_1 \
      spptrf_1 \
      spptrf_6 \
      spptri_1 \
      spptrs_1 \
      spteqr_1 \
      ssbtrd_1 \
      sspcon_1 \
      sspgst_1 \
      ssprfs_1 \
      ssptrd_1 \
      ssptrd_3 \
      ssptrf_1 \
      ssptri_1 \
      ssptrs_1 \
      sstebz_1 \
      sstebz_2 \
      sstebz_3 \
      sstebz_4 \
      sstein_1 \
      ssteqr_1 \
      ssteqr_3 \
      ssteqr_4 \
      ssterf_1 \
      ssterf_2 \
      ssterf_3 \
      ssterf_4 \
      ssterf_5 \
      ssycon_1 \
      ssygst_1 \
      ssyrfs_1 \
      ssytrd_1 \
      ssytrd_3 \
      ssytrf_1 \
      ssytri_1 \
      ssytrs_1 \
      stbcon_1 \
      stbrfs_1 \
      stbtrs_1 \
      stpcon_1 \
      stprfs_1 \
      stptri_1 \
      stptrs_1 \
      strcon_1 \
      strevc_1 \
      strevc_2 \
      strexc_1 \
      strrfs_1 \
      strsen_1 \
      strsna_1 \
      strsyl_1 \
      strtri_1 \
      strtrs_1 \
      zbdsqr_1 \
      zbdsqr_2 \
      zgbbrd_1 \
      zgbcon_1 \
      zgbrfs_1 \
      zgbtrf_1 \
      zgbtrs_1 \
      zgebak_1 \
      zgebal_1 \
      zgebrd_1 \
      zgebrd_3 \
      zgebrd_4 \
      zgebrd_5 \
      zgecon_1 \
      zgehrd_1 \
      zgehrd_2 \
      zgehrd_3 \
      zgelqf_1 \
      zgelqf_3 \
      zgeqpf_1 \
      zgeqrf_1 \
      zgeqrf_3 \
      zgerfs_1 \
      zgetrf_1 \
      zgetri_1 \
      zgetrs_1 \
      zhbtrd_1 \
      zhecon_1 \
      zhegst_1 \
      zherfs_1 \
      zhetrd_1 \
      zhetrd_2 \
      zhetrd_3 \
      zhetrf_1 \
      zhetri_1 \
      zhetrs_1 \
      zhpcon_1 \
      zhpgst_1 \
      zhprfs_1 \
      zhptrd_1 \
      zhptrd_2 \
      zhptrf_1 \
      zhptri_1 \
      zhptrs_1 \
      zhsein_1 \
      zhseqr_1 \
      zhseqr_2 \
      zhseqr_3 \
      zhseqr_4 \
      zpbcon_1 \
      zpbrfs_1 \
      zpbtrf_1 \
      zpbtrs_1 \
      zpocon_1 \
      zporfs_1 \
      zpotrf_1 \
      zpotri_1 \
      zpotrs_1 \
      zppcon_1 \
      zpprfs_1 \
      zpptrf_1 \
      zpptri_1 \
      zpptrs_1 \
      zpteqr_1 \
      zspcon_1 \
      zsprfs_1 \
      zsptrf_1 \
      zsptri_1 \
      zsptrs_1 \
      zstein_1 \
      zsteqr_1 \
      zsteqr_2 \
      zsycon_1 \
      zsyrfs_1 \
      zsytrf_1 \
      zsytri_1 \
      zsytrs_1 \
      ztbcon_1 \
      ztbrfs_1 \
      ztbtrs_1 \
      ztpcon_1 \
      ztprfs_1 \
      ztptri_1 \
      ztptrs_1 \
      ztrcon_1 \
      ztrevc_1 \
      ztrevc_2 \
      ztrexc_1 \
      ztrrfs_1 \
      ztrsen_1 \
      ztrsna_1 \
      ztrsyl_1 \
      ztrtri_1 \
      ztrtrs_1 \
      zungbr_1 \
      zungbr_2 \
      zungbr_3 \
      zungbr_4 \
      zunghr_1 \
      zunghr_2 \
      zunglq_1 \
      zunglq_2 \
      zungqr_1 \
      zungqr_2 \
      zungtr_1 \
      zungtr_2 \
      zunmbr_1 \
      zunmbr_2 \
      zunmhr_1 \
      zunmlq_1 \
      zunmqr_1 \
      zunmqr_2 \
      zunmtr_1 \
      zupgtr_1 \
      zupmtr_1

OBJ = $(TST:%=%.o)
EXE = $(TST:%=x%)
OUT = $(TST:%=%.out)

AUX = ../utils/compare_complex_doubles.o \
      ../utils/compare_complex_floats.o \
      ../utils/compare_doubles.o \
      ../utils/compare_floats.o

all: tst

tst: $(OUT)

%.out : x%
	./$< > $@

x% : %.o
	$(LINKER) $(LDFLAGS) $< -o $@ $(AUX) ../../$(LAPACKE) $(LIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) -I ../include -I ../../include -o $@ $<

clean:
	rm -f $(EXE) $(OBJ)
