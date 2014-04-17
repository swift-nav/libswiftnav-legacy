# Disable code coverage for this subdirectory
set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")

# Build libraries with PIC flag as they will be linked into shared objects
if (NOT CMAKE_CROSSCOMPILING)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(OBJ_UTILS utils/lapacke_cgb_nancheck.c utils/lapacke_cgb_trans.c
  utils/lapacke_cge_nancheck.c utils/lapacke_cge_trans.c
  utils/lapacke_cgg_nancheck.c utils/lapacke_cgg_trans.c
  utils/lapacke_cgt_nancheck.c utils/lapacke_chb_nancheck.c
  utils/lapacke_chb_trans.c utils/lapacke_che_nancheck.c
  utils/lapacke_che_trans.c utils/lapacke_chp_nancheck.c
  utils/lapacke_chp_trans.c utils/lapacke_chs_nancheck.c
  utils/lapacke_chs_trans.c utils/lapacke_c_nancheck.c
  utils/lapacke_cpb_nancheck.c utils/lapacke_cpb_trans.c
  utils/lapacke_cpf_nancheck.c utils/lapacke_cpf_trans.c
  utils/lapacke_cpo_nancheck.c utils/lapacke_cpo_trans.c
  utils/lapacke_cpp_nancheck.c utils/lapacke_cpp_trans.c
  utils/lapacke_cpt_nancheck.c utils/lapacke_csp_nancheck.c
  utils/lapacke_csp_trans.c utils/lapacke_cst_nancheck.c
  utils/lapacke_csy_nancheck.c utils/lapacke_csy_trans.c
  utils/lapacke_ctb_nancheck.c utils/lapacke_ctb_trans.c
  utils/lapacke_ctf_nancheck.c utils/lapacke_ctf_trans.c
  utils/lapacke_ctp_nancheck.c utils/lapacke_ctp_trans.c
  utils/lapacke_ctr_nancheck.c utils/lapacke_ctr_trans.c
  utils/lapacke_dgb_nancheck.c utils/lapacke_dgb_trans.c
  utils/lapacke_dge_nancheck.c utils/lapacke_dge_trans.c
  utils/lapacke_dgg_nancheck.c utils/lapacke_dgg_trans.c
  utils/lapacke_dgt_nancheck.c utils/lapacke_dhs_nancheck.c
  utils/lapacke_dhs_trans.c utils/lapacke_d_nancheck.c
  utils/lapacke_dpb_nancheck.c utils/lapacke_dpb_trans.c
  utils/lapacke_dpf_nancheck.c utils/lapacke_dpf_trans.c
  utils/lapacke_dpo_nancheck.c utils/lapacke_dpo_trans.c
  utils/lapacke_dpp_nancheck.c utils/lapacke_dpp_trans.c
  utils/lapacke_dpt_nancheck.c utils/lapacke_dsb_nancheck.c
  utils/lapacke_dsb_trans.c utils/lapacke_dsp_nancheck.c
  utils/lapacke_dsp_trans.c utils/lapacke_dst_nancheck.c
  utils/lapacke_dsy_nancheck.c utils/lapacke_dsy_trans.c
  utils/lapacke_dtb_nancheck.c utils/lapacke_dtb_trans.c
  utils/lapacke_dtf_nancheck.c utils/lapacke_dtf_trans.c
  utils/lapacke_dtp_nancheck.c utils/lapacke_dtp_trans.c
  utils/lapacke_dtr_nancheck.c utils/lapacke_dtr_trans.c utils/lapacke_lsame.c
  utils/lapacke_sgb_nancheck.c utils/lapacke_sgb_trans.c
  utils/lapacke_sge_nancheck.c utils/lapacke_sge_trans.c
  utils/lapacke_sgg_nancheck.c utils/lapacke_sgg_trans.c
  utils/lapacke_sgt_nancheck.c utils/lapacke_shs_nancheck.c
  utils/lapacke_shs_trans.c utils/lapacke_s_nancheck.c
  utils/lapacke_spb_nancheck.c utils/lapacke_spb_trans.c
  utils/lapacke_spf_nancheck.c utils/lapacke_spf_trans.c
  utils/lapacke_spo_nancheck.c utils/lapacke_spo_trans.c
  utils/lapacke_spp_nancheck.c utils/lapacke_spp_trans.c
  utils/lapacke_spt_nancheck.c utils/lapacke_ssb_nancheck.c
  utils/lapacke_ssb_trans.c utils/lapacke_ssp_nancheck.c
  utils/lapacke_ssp_trans.c utils/lapacke_sst_nancheck.c
  utils/lapacke_ssy_nancheck.c utils/lapacke_ssy_trans.c
  utils/lapacke_stb_nancheck.c utils/lapacke_stb_trans.c
  utils/lapacke_stf_nancheck.c utils/lapacke_stf_trans.c
  utils/lapacke_stp_nancheck.c utils/lapacke_stp_trans.c
  utils/lapacke_str_nancheck.c utils/lapacke_str_trans.c utils/lapacke_xerbla.c
  utils/lapacke_zgb_nancheck.c utils/lapacke_zgb_trans.c
  utils/lapacke_zge_nancheck.c utils/lapacke_zge_trans.c
  utils/lapacke_zgg_nancheck.c utils/lapacke_zgg_trans.c
  utils/lapacke_zgt_nancheck.c utils/lapacke_zhb_nancheck.c
  utils/lapacke_zhb_trans.c utils/lapacke_zhe_nancheck.c
  utils/lapacke_zhe_trans.c utils/lapacke_zhp_nancheck.c
  utils/lapacke_zhp_trans.c utils/lapacke_zhs_nancheck.c
  utils/lapacke_zhs_trans.c utils/lapacke_z_nancheck.c
  utils/lapacke_zpb_nancheck.c utils/lapacke_zpb_trans.c
  utils/lapacke_zpf_nancheck.c utils/lapacke_zpf_trans.c
  utils/lapacke_zpo_nancheck.c utils/lapacke_zpo_trans.c
  utils/lapacke_zpp_nancheck.c utils/lapacke_zpp_trans.c
  utils/lapacke_zpt_nancheck.c utils/lapacke_zsp_nancheck.c
  utils/lapacke_zsp_trans.c utils/lapacke_zst_nancheck.c
  utils/lapacke_zsy_nancheck.c utils/lapacke_zsy_trans.c
  utils/lapacke_ztb_nancheck.c utils/lapacke_ztb_trans.c
  utils/lapacke_ztf_nancheck.c utils/lapacke_ztf_trans.c
  utils/lapacke_ztp_nancheck.c utils/lapacke_ztp_trans.c
  utils/lapacke_ztr_nancheck.c utils/lapacke_ztr_trans.c
  utils/lapacke_make_complex_float.c utils/lapacke_make_complex_double.c)

set(OBJ src/lapacke_cbdsqr.c src/lapacke_cbdsqr_work.c src/lapacke_cgbbrd.c
  src/lapacke_cgbbrd_work.c src/lapacke_cgbcon.c src/lapacke_cgbcon_work.c
  src/lapacke_cgbequb.c src/lapacke_cgbequb_work.c src/lapacke_cgbequ.c
  src/lapacke_cgbequ_work.c src/lapacke_cgbrfs.c src/lapacke_cgbrfs_work.c
  src/lapacke_cgbrfsx.c src/lapacke_cgbrfsx_work.c src/lapacke_cgbsv.c
  src/lapacke_cgbsv_work.c src/lapacke_cgbsvx.c src/lapacke_cgbsvx_work.c
  src/lapacke_cgbsvxx.c src/lapacke_cgbsvxx_work.c src/lapacke_cgbtrf.c
  src/lapacke_cgbtrf_work.c src/lapacke_cgbtrs.c src/lapacke_cgbtrs_work.c
  src/lapacke_cgebak.c src/lapacke_cgebak_work.c src/lapacke_cgebal.c
  src/lapacke_cgebal_work.c src/lapacke_cgebrd.c src/lapacke_cgebrd_work.c
  src/lapacke_cgecon.c src/lapacke_cgecon_work.c src/lapacke_cgeequb.c
  src/lapacke_cgeequb_work.c src/lapacke_cgeequ.c src/lapacke_cgeequ_work.c
  src/lapacke_cgees.c src/lapacke_cgees_work.c src/lapacke_cgeesx.c
  src/lapacke_cgeesx_work.c src/lapacke_cgeev.c src/lapacke_cgeev_work.c
  src/lapacke_cgeevx.c src/lapacke_cgeevx_work.c src/lapacke_cgehrd.c
  src/lapacke_cgehrd_work.c src/lapacke_cgelqf.c src/lapacke_cgelqf_work.c
  src/lapacke_cgels.c src/lapacke_cgelsd.c src/lapacke_cgelsd_work.c
  src/lapacke_cgelss.c src/lapacke_cgelss_work.c src/lapacke_cgels_work.c
  src/lapacke_cgelsy.c src/lapacke_cgelsy_work.c src/lapacke_cgeqlf.c
  src/lapacke_cgeqlf_work.c src/lapacke_cgeqp3.c src/lapacke_cgeqp3_work.c
  src/lapacke_cgeqpf.c src/lapacke_cgeqpf_work.c src/lapacke_cgeqrf.c
  src/lapacke_cgeqrf_work.c src/lapacke_cgerfs.c src/lapacke_cgerfs_work.c
  src/lapacke_cgerfsx.c src/lapacke_cgerfsx_work.c src/lapacke_cgerqf.c
  src/lapacke_cgerqf_work.c src/lapacke_cgesdd.c src/lapacke_cgesdd_work.c
  src/lapacke_cgesv.c src/lapacke_cgesvd.c src/lapacke_cgesvd_work.c
  src/lapacke_cgesv_work.c src/lapacke_cgesvx.c src/lapacke_cgesvx_work.c
  src/lapacke_cgesvxx.c src/lapacke_cgesvxx_work.c src/lapacke_cgetrf.c
  src/lapacke_cgetrf_work.c src/lapacke_cgetri.c src/lapacke_cgetri_work.c
  src/lapacke_cgetrs.c src/lapacke_cgetrs_work.c src/lapacke_cggbak.c
  src/lapacke_cggbak_work.c src/lapacke_cggbal.c src/lapacke_cggbal_work.c
  src/lapacke_cgges.c src/lapacke_cgges_work.c src/lapacke_cggesx.c
  src/lapacke_cggesx_work.c src/lapacke_cggev.c src/lapacke_cggev_work.c
  src/lapacke_cggevx.c src/lapacke_cggevx_work.c src/lapacke_cggglm.c
  src/lapacke_cggglm_work.c src/lapacke_cgghrd.c src/lapacke_cgghrd_work.c
  src/lapacke_cgglse.c src/lapacke_cgglse_work.c src/lapacke_cggqrf.c
  src/lapacke_cggqrf_work.c src/lapacke_cggrqf.c src/lapacke_cggrqf_work.c
  src/lapacke_cggsvd.c src/lapacke_cggsvd_work.c src/lapacke_cggsvp.c
  src/lapacke_cggsvp_work.c src/lapacke_cgtcon.c src/lapacke_cgtcon_work.c
  src/lapacke_cgtrfs.c src/lapacke_cgtrfs_work.c src/lapacke_cgtsv.c
  src/lapacke_cgtsv_work.c src/lapacke_cgtsvx.c src/lapacke_cgtsvx_work.c
  src/lapacke_cgttrf.c src/lapacke_cgttrf_work.c src/lapacke_cgttrs.c
  src/lapacke_cgttrs_work.c src/lapacke_chbev.c src/lapacke_chbevd.c
  src/lapacke_chbevd_work.c src/lapacke_chbev_work.c src/lapacke_chbevx.c
  src/lapacke_chbevx_work.c src/lapacke_chbgst.c src/lapacke_chbgst_work.c
  src/lapacke_chbgv.c src/lapacke_chbgvd.c src/lapacke_chbgvd_work.c
  src/lapacke_chbgv_work.c src/lapacke_chbgvx.c src/lapacke_chbgvx_work.c
  src/lapacke_chbtrd.c src/lapacke_chbtrd_work.c src/lapacke_checon.c
  src/lapacke_checon_work.c src/lapacke_cheequb.c src/lapacke_cheequb_work.c
  src/lapacke_cheev.c src/lapacke_cheevd.c src/lapacke_cheevd_work.c
  src/lapacke_cheevr.c src/lapacke_cheevr_work.c src/lapacke_cheev_work.c
  src/lapacke_cheevx.c src/lapacke_cheevx_work.c src/lapacke_chegst.c
  src/lapacke_chegst_work.c src/lapacke_chegv.c src/lapacke_chegvd.c
  src/lapacke_chegvd_work.c src/lapacke_chegv_work.c src/lapacke_chegvx.c
  src/lapacke_chegvx_work.c src/lapacke_cherfs.c src/lapacke_cherfs_work.c
  src/lapacke_cherfsx.c src/lapacke_cherfsx_work.c src/lapacke_chesv.c
  src/lapacke_chesv_work.c src/lapacke_chesvx.c src/lapacke_chesvx_work.c
  src/lapacke_chesvxx.c src/lapacke_chesvxx_work.c src/lapacke_chetrd.c
  src/lapacke_chetrd_work.c src/lapacke_chetrf.c src/lapacke_chetrf_work.c
  src/lapacke_chetri.c src/lapacke_chetri_work.c src/lapacke_chetrs.c
  src/lapacke_chetrs_work.c src/lapacke_chfrk.c src/lapacke_chfrk_work.c
  src/lapacke_chgeqz.c src/lapacke_chgeqz_work.c src/lapacke_chpcon.c
  src/lapacke_chpcon_work.c src/lapacke_chpev.c src/lapacke_chpevd.c
  src/lapacke_chpevd_work.c src/lapacke_chpev_work.c src/lapacke_chpevx.c
  src/lapacke_chpevx_work.c src/lapacke_chpgst.c src/lapacke_chpgst_work.c
  src/lapacke_chpgv.c src/lapacke_chpgvd.c src/lapacke_chpgvd_work.c
  src/lapacke_chpgv_work.c src/lapacke_chpgvx.c src/lapacke_chpgvx_work.c
  src/lapacke_chprfs.c src/lapacke_chprfs_work.c src/lapacke_chpsv.c
  src/lapacke_chpsv_work.c src/lapacke_chpsvx.c src/lapacke_chpsvx_work.c
  src/lapacke_chptrd.c src/lapacke_chptrd_work.c src/lapacke_chptrf.c
  src/lapacke_chptrf_work.c src/lapacke_chptri.c src/lapacke_chptri_work.c
  src/lapacke_chptrs.c src/lapacke_chptrs_work.c src/lapacke_chsein.c
  src/lapacke_chsein_work.c src/lapacke_chseqr.c src/lapacke_chseqr_work.c
  src/lapacke_cpbcon.c src/lapacke_cpbcon_work.c src/lapacke_cpbequ.c
  src/lapacke_cpbequ_work.c src/lapacke_cpbrfs.c src/lapacke_cpbrfs_work.c
  src/lapacke_cpbstf.c src/lapacke_cpbstf_work.c src/lapacke_cpbsv.c
  src/lapacke_cpbsv_work.c src/lapacke_cpbsvx.c src/lapacke_cpbsvx_work.c
  src/lapacke_cpbtrf.c src/lapacke_cpbtrf_work.c src/lapacke_cpbtrs.c
  src/lapacke_cpbtrs_work.c src/lapacke_cpftrf.c src/lapacke_cpftrf_work.c
  src/lapacke_cpftri.c src/lapacke_cpftri_work.c src/lapacke_cpftrs.c
  src/lapacke_cpftrs_work.c src/lapacke_cpocon.c src/lapacke_cpocon_work.c
  src/lapacke_cpoequb.c src/lapacke_cpoequb_work.c src/lapacke_cpoequ.c
  src/lapacke_cpoequ_work.c src/lapacke_cporfs.c src/lapacke_cporfs_work.c
  src/lapacke_cporfsx.c src/lapacke_cporfsx_work.c src/lapacke_cposv.c
  src/lapacke_cposv_work.c src/lapacke_cposvx.c src/lapacke_cposvx_work.c
  src/lapacke_cposvxx.c src/lapacke_cposvxx_work.c src/lapacke_cpotrf.c
  src/lapacke_cpotrf_work.c src/lapacke_cpotri.c src/lapacke_cpotri_work.c
  src/lapacke_cpotrs.c src/lapacke_cpotrs_work.c src/lapacke_cppcon.c
  src/lapacke_cppcon_work.c src/lapacke_cppequ.c src/lapacke_cppequ_work.c
  src/lapacke_cpprfs.c src/lapacke_cpprfs_work.c src/lapacke_cppsv.c
  src/lapacke_cppsv_work.c src/lapacke_cppsvx.c src/lapacke_cppsvx_work.c
  src/lapacke_cpptrf.c src/lapacke_cpptrf_work.c src/lapacke_cpptri.c
  src/lapacke_cpptri_work.c src/lapacke_cpptrs.c src/lapacke_cpptrs_work.c
  src/lapacke_cpstrf.c src/lapacke_cpstrf_work.c src/lapacke_cptcon.c
  src/lapacke_cptcon_work.c src/lapacke_cpteqr.c src/lapacke_cpteqr_work.c
  src/lapacke_cptrfs.c src/lapacke_cptrfs_work.c src/lapacke_cptsv.c
  src/lapacke_cptsv_work.c src/lapacke_cptsvx.c src/lapacke_cptsvx_work.c
  src/lapacke_cpttrf.c src/lapacke_cpttrf_work.c src/lapacke_cpttrs.c
  src/lapacke_cpttrs_work.c src/lapacke_cspcon.c src/lapacke_cspcon_work.c
  src/lapacke_csprfs.c src/lapacke_csprfs_work.c src/lapacke_cspsv.c
  src/lapacke_cspsv_work.c src/lapacke_cspsvx.c src/lapacke_cspsvx_work.c
  src/lapacke_csptrf.c src/lapacke_csptrf_work.c src/lapacke_csptri.c
  src/lapacke_csptri_work.c src/lapacke_csptrs.c src/lapacke_csptrs_work.c
  src/lapacke_cstedc.c src/lapacke_cstedc_work.c src/lapacke_cstegr.c
  src/lapacke_cstegr_work.c src/lapacke_cstein.c src/lapacke_cstein_work.c
  src/lapacke_cstemr.c src/lapacke_cstemr_work.c src/lapacke_csteqr.c
  src/lapacke_csteqr_work.c src/lapacke_csycon.c src/lapacke_csycon_work.c
  src/lapacke_csyequb.c src/lapacke_csyequb_work.c src/lapacke_csyrfs.c
  src/lapacke_csyrfs_work.c src/lapacke_csyrfsx.c src/lapacke_csyrfsx_work.c
  src/lapacke_csysv.c src/lapacke_csysv_work.c src/lapacke_csysvx.c
  src/lapacke_csysvx_work.c src/lapacke_csysvxx.c src/lapacke_csysvxx_work.c
  src/lapacke_csytrf.c src/lapacke_csytrf_work.c src/lapacke_csytri.c
  src/lapacke_csytri_work.c src/lapacke_csytrs.c src/lapacke_csytrs_work.c
  src/lapacke_ctbcon.c src/lapacke_ctbcon_work.c src/lapacke_ctbrfs.c
  src/lapacke_ctbrfs_work.c src/lapacke_ctbtrs.c src/lapacke_ctbtrs_work.c
  src/lapacke_ctfsm.c src/lapacke_ctfsm_work.c src/lapacke_ctftri.c
  src/lapacke_ctftri_work.c src/lapacke_ctfttp.c src/lapacke_ctfttp_work.c
  src/lapacke_ctfttr.c src/lapacke_ctfttr_work.c src/lapacke_ctgevc.c
  src/lapacke_ctgevc_work.c src/lapacke_ctgexc.c src/lapacke_ctgexc_work.c
  src/lapacke_ctgsen.c src/lapacke_ctgsen_work.c src/lapacke_ctgsja.c
  src/lapacke_ctgsja_work.c src/lapacke_ctgsna.c src/lapacke_ctgsna_work.c
  src/lapacke_ctgsyl.c src/lapacke_ctgsyl_work.c src/lapacke_ctpcon.c
  src/lapacke_ctpcon_work.c src/lapacke_ctprfs.c src/lapacke_ctprfs_work.c
  src/lapacke_ctptri.c src/lapacke_ctptri_work.c src/lapacke_ctptrs.c
  src/lapacke_ctptrs_work.c src/lapacke_ctpttf.c src/lapacke_ctpttf_work.c
  src/lapacke_ctpttr.c src/lapacke_ctpttr_work.c src/lapacke_ctrcon.c
  src/lapacke_ctrcon_work.c src/lapacke_ctrevc.c src/lapacke_ctrevc_work.c
  src/lapacke_ctrexc.c src/lapacke_ctrexc_work.c src/lapacke_ctrrfs.c
  src/lapacke_ctrrfs_work.c src/lapacke_ctrsen.c src/lapacke_ctrsen_work.c
  src/lapacke_ctrsna.c src/lapacke_ctrsna_work.c src/lapacke_ctrsyl.c
  src/lapacke_ctrsyl_work.c src/lapacke_ctrtri.c src/lapacke_ctrtri_work.c
  src/lapacke_ctrtrs.c src/lapacke_ctrtrs_work.c src/lapacke_ctrttf.c
  src/lapacke_ctrttf_work.c src/lapacke_ctrttp.c src/lapacke_ctrttp_work.c
  src/lapacke_ctzrzf.c src/lapacke_ctzrzf_work.c src/lapacke_cungbr.c
  src/lapacke_cungbr_work.c src/lapacke_cunghr.c src/lapacke_cunghr_work.c
  src/lapacke_cunglq.c src/lapacke_cunglq_work.c src/lapacke_cungql.c
  src/lapacke_cungql_work.c src/lapacke_cungqr.c src/lapacke_cungqr_work.c
  src/lapacke_cungrq.c src/lapacke_cungrq_work.c src/lapacke_cungtr.c
  src/lapacke_cungtr_work.c src/lapacke_cunmbr.c src/lapacke_cunmbr_work.c
  src/lapacke_cunmhr.c src/lapacke_cunmhr_work.c src/lapacke_cunmlq.c
  src/lapacke_cunmlq_work.c src/lapacke_cunmql.c src/lapacke_cunmql_work.c
  src/lapacke_cunmqr.c src/lapacke_cunmqr_work.c src/lapacke_cunmrq.c
  src/lapacke_cunmrq_work.c src/lapacke_cunmrz.c src/lapacke_cunmrz_work.c
  src/lapacke_cunmtr.c src/lapacke_cunmtr_work.c src/lapacke_cupgtr.c
  src/lapacke_cupgtr_work.c src/lapacke_cupmtr.c src/lapacke_cupmtr_work.c
  src/lapacke_dbdsdc.c src/lapacke_dbdsdc_work.c src/lapacke_dbdsqr.c
  src/lapacke_dbdsqr_work.c src/lapacke_ddisna.c src/lapacke_ddisna_work.c
  src/lapacke_dgbbrd.c src/lapacke_dgbbrd_work.c src/lapacke_dgbcon.c
  src/lapacke_dgbcon_work.c src/lapacke_dgbequb.c src/lapacke_dgbequb_work.c
  src/lapacke_dgbequ.c src/lapacke_dgbequ_work.c src/lapacke_dgbrfs.c
  src/lapacke_dgbrfs_work.c src/lapacke_dgbrfsx.c src/lapacke_dgbrfsx_work.c
  src/lapacke_dgbsv.c src/lapacke_dgbsv_work.c src/lapacke_dgbsvx.c
  src/lapacke_dgbsvx_work.c src/lapacke_dgbsvxx.c src/lapacke_dgbsvxx_work.c
  src/lapacke_dgbtrf.c src/lapacke_dgbtrf_work.c src/lapacke_dgbtrs.c
  src/lapacke_dgbtrs_work.c src/lapacke_dgebak.c src/lapacke_dgebak_work.c
  src/lapacke_dgebal.c src/lapacke_dgebal_work.c src/lapacke_dgebrd.c
  src/lapacke_dgebrd_work.c src/lapacke_dgecon.c src/lapacke_dgecon_work.c
  src/lapacke_dgeequb.c src/lapacke_dgeequb_work.c src/lapacke_dgeequ.c
  src/lapacke_dgeequ_work.c src/lapacke_dgees.c src/lapacke_dgees_work.c
  src/lapacke_dgeesx.c src/lapacke_dgeesx_work.c src/lapacke_dgeev.c
  src/lapacke_dgeev_work.c src/lapacke_dgeevx.c src/lapacke_dgeevx_work.c
  src/lapacke_dgehrd.c src/lapacke_dgehrd_work.c src/lapacke_dgejsv.c
  src/lapacke_dgejsv_work.c src/lapacke_dgelqf.c src/lapacke_dgelqf_work.c
  src/lapacke_dgels.c src/lapacke_dgelsd.c src/lapacke_dgelsd_work.c
  src/lapacke_dgelss.c src/lapacke_dgelss_work.c src/lapacke_dgels_work.c
  src/lapacke_dgelsy.c src/lapacke_dgelsy_work.c src/lapacke_dgeqlf.c
  src/lapacke_dgeqlf_work.c src/lapacke_dgeqp3.c src/lapacke_dgeqp3_work.c
  src/lapacke_dgeqpf.c src/lapacke_dgeqpf_work.c src/lapacke_dgeqrf.c
  src/lapacke_dgeqrf_work.c src/lapacke_dgerfs.c src/lapacke_dgerfs_work.c
  src/lapacke_dgerfsx.c src/lapacke_dgerfsx_work.c src/lapacke_dgerqf.c
  src/lapacke_dgerqf_work.c src/lapacke_dgesdd.c src/lapacke_dgesdd_work.c
  src/lapacke_dgesv.c src/lapacke_dgesvd.c src/lapacke_dgesvd_work.c
  src/lapacke_dgesvj.c src/lapacke_dgesvj_work.c src/lapacke_dgesv_work.c
  src/lapacke_dgesvx.c src/lapacke_dgesvx_work.c src/lapacke_dgesvxx.c
  src/lapacke_dgesvxx_work.c src/lapacke_dgetrf.c src/lapacke_dgetrf_work.c
  src/lapacke_dgetri.c src/lapacke_dgetri_work.c src/lapacke_dgetrs.c
  src/lapacke_dgetrs_work.c src/lapacke_dggbak.c src/lapacke_dggbak_work.c
  src/lapacke_dggbal.c src/lapacke_dggbal_work.c src/lapacke_dgges.c
  src/lapacke_dgges_work.c src/lapacke_dggesx.c src/lapacke_dggesx_work.c
  src/lapacke_dggev.c src/lapacke_dggev_work.c src/lapacke_dggevx.c
  src/lapacke_dggevx_work.c src/lapacke_dggglm.c src/lapacke_dggglm_work.c
  src/lapacke_dgghrd.c src/lapacke_dgghrd_work.c src/lapacke_dgglse.c
  src/lapacke_dgglse_work.c src/lapacke_dggqrf.c src/lapacke_dggqrf_work.c
  src/lapacke_dggrqf.c src/lapacke_dggrqf_work.c src/lapacke_dggsvd.c
  src/lapacke_dggsvd_work.c src/lapacke_dggsvp.c src/lapacke_dggsvp_work.c
  src/lapacke_dgtcon.c src/lapacke_dgtcon_work.c src/lapacke_dgtrfs.c
  src/lapacke_dgtrfs_work.c src/lapacke_dgtsv.c src/lapacke_dgtsv_work.c
  src/lapacke_dgtsvx.c src/lapacke_dgtsvx_work.c src/lapacke_dgttrf.c
  src/lapacke_dgttrf_work.c src/lapacke_dgttrs.c src/lapacke_dgttrs_work.c
  src/lapacke_dhgeqz.c src/lapacke_dhgeqz_work.c src/lapacke_dhsein.c
  src/lapacke_dhsein_work.c src/lapacke_dhseqr.c src/lapacke_dhseqr_work.c
  src/lapacke_dopgtr.c src/lapacke_dopgtr_work.c src/lapacke_dopmtr.c
  src/lapacke_dopmtr_work.c src/lapacke_dorgbr.c src/lapacke_dorgbr_work.c
  src/lapacke_dorghr.c src/lapacke_dorghr_work.c src/lapacke_dorglq.c
  src/lapacke_dorglq_work.c src/lapacke_dorgql.c src/lapacke_dorgql_work.c
  src/lapacke_dorgqr.c src/lapacke_dorgqr_work.c src/lapacke_dorgrq.c
  src/lapacke_dorgrq_work.c src/lapacke_dorgtr.c src/lapacke_dorgtr_work.c
  src/lapacke_dormbr.c src/lapacke_dormbr_work.c src/lapacke_dormhr.c
  src/lapacke_dormhr_work.c src/lapacke_dormlq.c src/lapacke_dormlq_work.c
  src/lapacke_dormql.c src/lapacke_dormql_work.c src/lapacke_dormqr.c
  src/lapacke_dormqr_work.c src/lapacke_dormrq.c src/lapacke_dormrq_work.c
  src/lapacke_dormrz.c src/lapacke_dormrz_work.c src/lapacke_dormtr.c
  src/lapacke_dormtr_work.c src/lapacke_dpbcon.c src/lapacke_dpbcon_work.c
  src/lapacke_dpbequ.c src/lapacke_dpbequ_work.c src/lapacke_dpbrfs.c
  src/lapacke_dpbrfs_work.c src/lapacke_dpbstf.c src/lapacke_dpbstf_work.c
  src/lapacke_dpbsv.c src/lapacke_dpbsv_work.c src/lapacke_dpbsvx.c
  src/lapacke_dpbsvx_work.c src/lapacke_dpbtrf.c src/lapacke_dpbtrf_work.c
  src/lapacke_dpbtrs.c src/lapacke_dpbtrs_work.c src/lapacke_dpftrf.c
  src/lapacke_dpftrf_work.c src/lapacke_dpftri.c src/lapacke_dpftri_work.c
  src/lapacke_dpftrs.c src/lapacke_dpftrs_work.c src/lapacke_dpocon.c
  src/lapacke_dpocon_work.c src/lapacke_dpoequb.c src/lapacke_dpoequb_work.c
  src/lapacke_dpoequ.c src/lapacke_dpoequ_work.c src/lapacke_dporfs.c
  src/lapacke_dporfs_work.c src/lapacke_dporfsx.c src/lapacke_dporfsx_work.c
  src/lapacke_dposv.c src/lapacke_dposv_work.c src/lapacke_dposvx.c
  src/lapacke_dposvx_work.c src/lapacke_dposvxx.c src/lapacke_dposvxx_work.c
  src/lapacke_dpotrf.c src/lapacke_dpotrf_work.c src/lapacke_dpotri.c
  src/lapacke_dpotri_work.c src/lapacke_dpotrs.c src/lapacke_dpotrs_work.c
  src/lapacke_dppcon.c src/lapacke_dppcon_work.c src/lapacke_dppequ.c
  src/lapacke_dppequ_work.c src/lapacke_dpprfs.c src/lapacke_dpprfs_work.c
  src/lapacke_dppsv.c src/lapacke_dppsv_work.c src/lapacke_dppsvx.c
  src/lapacke_dppsvx_work.c src/lapacke_dpptrf.c src/lapacke_dpptrf_work.c
  src/lapacke_dpptri.c src/lapacke_dpptri_work.c src/lapacke_dpptrs.c
  src/lapacke_dpptrs_work.c src/lapacke_dpstrf.c src/lapacke_dpstrf_work.c
  src/lapacke_dptcon.c src/lapacke_dptcon_work.c src/lapacke_dpteqr.c
  src/lapacke_dpteqr_work.c src/lapacke_dptrfs.c src/lapacke_dptrfs_work.c
  src/lapacke_dptsv.c src/lapacke_dptsv_work.c src/lapacke_dptsvx.c
  src/lapacke_dptsvx_work.c src/lapacke_dpttrf.c src/lapacke_dpttrf_work.c
  src/lapacke_dpttrs.c src/lapacke_dpttrs_work.c src/lapacke_dsbev.c
  src/lapacke_dsbevd.c src/lapacke_dsbevd_work.c src/lapacke_dsbev_work.c
  src/lapacke_dsbevx.c src/lapacke_dsbevx_work.c src/lapacke_dsbgst.c
  src/lapacke_dsbgst_work.c src/lapacke_dsbgv.c src/lapacke_dsbgvd.c
  src/lapacke_dsbgvd_work.c src/lapacke_dsbgv_work.c src/lapacke_dsbgvx.c
  src/lapacke_dsbgvx_work.c src/lapacke_dsbtrd.c src/lapacke_dsbtrd_work.c
  src/lapacke_dsfrk.c src/lapacke_dsfrk_work.c src/lapacke_dsgesv.c
  src/lapacke_dsgesv_work.c src/lapacke_dspcon.c src/lapacke_dspcon_work.c
  src/lapacke_dspev.c src/lapacke_dspevd.c src/lapacke_dspevd_work.c
  src/lapacke_dspev_work.c src/lapacke_dspevx.c src/lapacke_dspevx_work.c
  src/lapacke_dspgst.c src/lapacke_dspgst_work.c src/lapacke_dspgv.c
  src/lapacke_dspgvd.c src/lapacke_dspgvd_work.c src/lapacke_dspgv_work.c
  src/lapacke_dspgvx.c src/lapacke_dspgvx_work.c src/lapacke_dsposv.c
  src/lapacke_dsposv_work.c src/lapacke_dsprfs.c src/lapacke_dsprfs_work.c
  src/lapacke_dspsv.c src/lapacke_dspsv_work.c src/lapacke_dspsvx.c
  src/lapacke_dspsvx_work.c src/lapacke_dsptrd.c src/lapacke_dsptrd_work.c
  src/lapacke_dsptrf.c src/lapacke_dsptrf_work.c src/lapacke_dsptri.c
  src/lapacke_dsptri_work.c src/lapacke_dsptrs.c src/lapacke_dsptrs_work.c
  src/lapacke_dstebz.c src/lapacke_dstebz_work.c src/lapacke_dstedc.c
  src/lapacke_dstedc_work.c src/lapacke_dstegr.c src/lapacke_dstegr_work.c
  src/lapacke_dstein.c src/lapacke_dstein_work.c src/lapacke_dstemr.c
  src/lapacke_dstemr_work.c src/lapacke_dsteqr.c src/lapacke_dsteqr_work.c
  src/lapacke_dsterf.c src/lapacke_dsterf_work.c src/lapacke_dstev.c
  src/lapacke_dstevd.c src/lapacke_dstevd_work.c src/lapacke_dstevr.c
  src/lapacke_dstevr_work.c src/lapacke_dstev_work.c src/lapacke_dstevx.c
  src/lapacke_dstevx_work.c src/lapacke_dsycon.c src/lapacke_dsycon_work.c
  src/lapacke_dsyequb.c src/lapacke_dsyequb_work.c src/lapacke_dsyev.c
  src/lapacke_dsyevd.c src/lapacke_dsyevd_work.c src/lapacke_dsyevr.c
  src/lapacke_dsyevr_work.c src/lapacke_dsyev_work.c src/lapacke_dsyevx.c
  src/lapacke_dsyevx_work.c src/lapacke_dsygst.c src/lapacke_dsygst_work.c
  src/lapacke_dsygv.c src/lapacke_dsygvd.c src/lapacke_dsygvd_work.c
  src/lapacke_dsygv_work.c src/lapacke_dsygvx.c src/lapacke_dsygvx_work.c
  src/lapacke_dsyrfs.c src/lapacke_dsyrfs_work.c src/lapacke_dsyrfsx.c
  src/lapacke_dsyrfsx_work.c src/lapacke_dsysv.c src/lapacke_dsysv_work.c
  src/lapacke_dsysvx.c src/lapacke_dsysvx_work.c src/lapacke_dsysvxx.c
  src/lapacke_dsysvxx_work.c src/lapacke_dsytrd.c src/lapacke_dsytrd_work.c
  src/lapacke_dsytrf.c src/lapacke_dsytrf_work.c src/lapacke_dsytri.c
  src/lapacke_dsytri_work.c src/lapacke_dsytrs.c src/lapacke_dsytrs_work.c
  src/lapacke_dtbcon.c src/lapacke_dtbcon_work.c src/lapacke_dtbrfs.c
  src/lapacke_dtbrfs_work.c src/lapacke_dtbtrs.c src/lapacke_dtbtrs_work.c
  src/lapacke_dtfsm.c src/lapacke_dtfsm_work.c src/lapacke_dtftri.c
  src/lapacke_dtftri_work.c src/lapacke_dtfttp.c src/lapacke_dtfttp_work.c
  src/lapacke_dtfttr.c src/lapacke_dtfttr_work.c src/lapacke_dtgevc.c
  src/lapacke_dtgevc_work.c src/lapacke_dtgexc.c src/lapacke_dtgexc_work.c
  src/lapacke_dtgsen.c src/lapacke_dtgsen_work.c src/lapacke_dtgsja.c
  src/lapacke_dtgsja_work.c src/lapacke_dtgsna.c src/lapacke_dtgsna_work.c
  src/lapacke_dtgsyl.c src/lapacke_dtgsyl_work.c src/lapacke_dtpcon.c
  src/lapacke_dtpcon_work.c src/lapacke_dtprfs.c src/lapacke_dtprfs_work.c
  src/lapacke_dtptri.c src/lapacke_dtptri_work.c src/lapacke_dtptrs.c
  src/lapacke_dtptrs_work.c src/lapacke_dtpttf.c src/lapacke_dtpttf_work.c
  src/lapacke_dtpttr.c src/lapacke_dtpttr_work.c src/lapacke_dtrcon.c
  src/lapacke_dtrcon_work.c src/lapacke_dtrevc.c src/lapacke_dtrevc_work.c
  src/lapacke_dtrexc.c src/lapacke_dtrexc_work.c src/lapacke_dtrrfs.c
  src/lapacke_dtrrfs_work.c src/lapacke_dtrsen.c src/lapacke_dtrsen_work.c
  src/lapacke_dtrsna.c src/lapacke_dtrsna_work.c src/lapacke_dtrsyl.c
  src/lapacke_dtrsyl_work.c src/lapacke_dtrtri.c src/lapacke_dtrtri_work.c
  src/lapacke_dtrtrs.c src/lapacke_dtrtrs_work.c src/lapacke_dtrttf.c
  src/lapacke_dtrttf_work.c src/lapacke_dtrttp.c src/lapacke_dtrttp_work.c
  src/lapacke_dtzrzf.c src/lapacke_dtzrzf_work.c src/lapacke_sbdsdc.c
  src/lapacke_sbdsdc_work.c src/lapacke_sbdsqr.c src/lapacke_sbdsqr_work.c
  src/lapacke_sdisna.c src/lapacke_sdisna_work.c src/lapacke_sgbbrd.c
  src/lapacke_sgbbrd_work.c src/lapacke_sgbcon.c src/lapacke_sgbcon_work.c
  src/lapacke_sgbequb.c src/lapacke_sgbequb_work.c src/lapacke_sgbequ.c
  src/lapacke_sgbequ_work.c src/lapacke_sgbrfs.c src/lapacke_sgbrfs_work.c
  src/lapacke_sgbrfsx.c src/lapacke_sgbrfsx_work.c src/lapacke_sgbsv.c
  src/lapacke_sgbsv_work.c src/lapacke_sgbsvx.c src/lapacke_sgbsvx_work.c
  src/lapacke_sgbsvxx.c src/lapacke_sgbsvxx_work.c src/lapacke_sgbtrf.c
  src/lapacke_sgbtrf_work.c src/lapacke_sgbtrs.c src/lapacke_sgbtrs_work.c
  src/lapacke_sgebak.c src/lapacke_sgebak_work.c src/lapacke_sgebal.c
  src/lapacke_sgebal_work.c src/lapacke_sgebrd.c src/lapacke_sgebrd_work.c
  src/lapacke_sgecon.c src/lapacke_sgecon_work.c src/lapacke_sgeequb.c
  src/lapacke_sgeequb_work.c src/lapacke_sgeequ.c src/lapacke_sgeequ_work.c
  src/lapacke_sgees.c src/lapacke_sgees_work.c src/lapacke_sgeesx.c
  src/lapacke_sgeesx_work.c src/lapacke_sgeev.c src/lapacke_sgeev_work.c
  src/lapacke_sgeevx.c src/lapacke_sgeevx_work.c src/lapacke_sgehrd.c
  src/lapacke_sgehrd_work.c src/lapacke_sgejsv.c src/lapacke_sgejsv_work.c
  src/lapacke_sgelqf.c src/lapacke_sgelqf_work.c src/lapacke_sgels.c
  src/lapacke_sgelsd.c src/lapacke_sgelsd_work.c src/lapacke_sgelss.c
  src/lapacke_sgelss_work.c src/lapacke_sgels_work.c src/lapacke_sgelsy.c
  src/lapacke_sgelsy_work.c src/lapacke_sgeqlf.c src/lapacke_sgeqlf_work.c
  src/lapacke_sgeqp3.c src/lapacke_sgeqp3_work.c src/lapacke_sgeqpf.c
  src/lapacke_sgeqpf_work.c src/lapacke_sgeqrf.c src/lapacke_sgeqrf_work.c
  src/lapacke_sgerfs.c src/lapacke_sgerfs_work.c src/lapacke_sgerfsx.c
  src/lapacke_sgerfsx_work.c src/lapacke_sgerqf.c src/lapacke_sgerqf_work.c
  src/lapacke_sgesdd.c src/lapacke_sgesdd_work.c src/lapacke_sgesv.c
  src/lapacke_sgesvd.c src/lapacke_sgesvd_work.c src/lapacke_sgesvj.c
  src/lapacke_sgesvj_work.c src/lapacke_sgesv_work.c src/lapacke_sgesvx.c
  src/lapacke_sgesvx_work.c src/lapacke_sgesvxx.c src/lapacke_sgesvxx_work.c
  src/lapacke_sgetrf.c src/lapacke_sgetrf_work.c src/lapacke_sgetri.c
  src/lapacke_sgetri_work.c src/lapacke_sgetrs.c src/lapacke_sgetrs_work.c
  src/lapacke_sggbak.c src/lapacke_sggbak_work.c src/lapacke_sggbal.c
  src/lapacke_sggbal_work.c src/lapacke_sgges.c src/lapacke_sgges_work.c
  src/lapacke_sggesx.c src/lapacke_sggesx_work.c src/lapacke_sggev.c
  src/lapacke_sggev_work.c src/lapacke_sggevx.c src/lapacke_sggevx_work.c
  src/lapacke_sggglm.c src/lapacke_sggglm_work.c src/lapacke_sgghrd.c
  src/lapacke_sgghrd_work.c src/lapacke_sgglse.c src/lapacke_sgglse_work.c
  src/lapacke_sggqrf.c src/lapacke_sggqrf_work.c src/lapacke_sggrqf.c
  src/lapacke_sggrqf_work.c src/lapacke_sggsvd.c src/lapacke_sggsvd_work.c
  src/lapacke_sggsvp.c src/lapacke_sggsvp_work.c src/lapacke_sgtcon.c
  src/lapacke_sgtcon_work.c src/lapacke_sgtrfs.c src/lapacke_sgtrfs_work.c
  src/lapacke_sgtsv.c src/lapacke_sgtsv_work.c src/lapacke_sgtsvx.c
  src/lapacke_sgtsvx_work.c src/lapacke_sgttrf.c src/lapacke_sgttrf_work.c
  src/lapacke_sgttrs.c src/lapacke_sgttrs_work.c src/lapacke_shgeqz.c
  src/lapacke_shgeqz_work.c src/lapacke_shsein.c src/lapacke_shsein_work.c
  src/lapacke_shseqr.c src/lapacke_shseqr_work.c src/lapacke_sopgtr.c
  src/lapacke_sopgtr_work.c src/lapacke_sopmtr.c src/lapacke_sopmtr_work.c
  src/lapacke_sorgbr.c src/lapacke_sorgbr_work.c src/lapacke_sorghr.c
  src/lapacke_sorghr_work.c src/lapacke_sorglq.c src/lapacke_sorglq_work.c
  src/lapacke_sorgql.c src/lapacke_sorgql_work.c src/lapacke_sorgqr.c
  src/lapacke_sorgqr_work.c src/lapacke_sorgrq.c src/lapacke_sorgrq_work.c
  src/lapacke_sorgtr.c src/lapacke_sorgtr_work.c src/lapacke_sormbr.c
  src/lapacke_sormbr_work.c src/lapacke_sormhr.c src/lapacke_sormhr_work.c
  src/lapacke_sormlq.c src/lapacke_sormlq_work.c src/lapacke_sormql.c
  src/lapacke_sormql_work.c src/lapacke_sormqr.c src/lapacke_sormqr_work.c
  src/lapacke_sormrq.c src/lapacke_sormrq_work.c src/lapacke_sormrz.c
  src/lapacke_sormrz_work.c src/lapacke_sormtr.c src/lapacke_sormtr_work.c
  src/lapacke_spbcon.c src/lapacke_spbcon_work.c src/lapacke_spbequ.c
  src/lapacke_spbequ_work.c src/lapacke_spbrfs.c src/lapacke_spbrfs_work.c
  src/lapacke_spbstf.c src/lapacke_spbstf_work.c src/lapacke_spbsv.c
  src/lapacke_spbsv_work.c src/lapacke_spbsvx.c src/lapacke_spbsvx_work.c
  src/lapacke_spbtrf.c src/lapacke_spbtrf_work.c src/lapacke_spbtrs.c
  src/lapacke_spbtrs_work.c src/lapacke_spftrf.c src/lapacke_spftrf_work.c
  src/lapacke_spftri.c src/lapacke_spftri_work.c src/lapacke_spftrs.c
  src/lapacke_spftrs_work.c src/lapacke_spocon.c src/lapacke_spocon_work.c
  src/lapacke_spoequb.c src/lapacke_spoequb_work.c src/lapacke_spoequ.c
  src/lapacke_spoequ_work.c src/lapacke_sporfs.c src/lapacke_sporfs_work.c
  src/lapacke_sporfsx.c src/lapacke_sporfsx_work.c src/lapacke_sposv.c
  src/lapacke_sposv_work.c src/lapacke_sposvx.c src/lapacke_sposvx_work.c
  src/lapacke_sposvxx.c src/lapacke_sposvxx_work.c src/lapacke_spotrf.c
  src/lapacke_spotrf_work.c src/lapacke_spotri.c src/lapacke_spotri_work.c
  src/lapacke_spotrs.c src/lapacke_spotrs_work.c src/lapacke_sppcon.c
  src/lapacke_sppcon_work.c src/lapacke_sppequ.c src/lapacke_sppequ_work.c
  src/lapacke_spprfs.c src/lapacke_spprfs_work.c src/lapacke_sppsv.c
  src/lapacke_sppsv_work.c src/lapacke_sppsvx.c src/lapacke_sppsvx_work.c
  src/lapacke_spptrf.c src/lapacke_spptrf_work.c src/lapacke_spptri.c
  src/lapacke_spptri_work.c src/lapacke_spptrs.c src/lapacke_spptrs_work.c
  src/lapacke_spstrf.c src/lapacke_spstrf_work.c src/lapacke_sptcon.c
  src/lapacke_sptcon_work.c src/lapacke_spteqr.c src/lapacke_spteqr_work.c
  src/lapacke_sptrfs.c src/lapacke_sptrfs_work.c src/lapacke_sptsv.c
  src/lapacke_sptsv_work.c src/lapacke_sptsvx.c src/lapacke_sptsvx_work.c
  src/lapacke_spttrf.c src/lapacke_spttrf_work.c src/lapacke_spttrs.c
  src/lapacke_spttrs_work.c src/lapacke_ssbev.c src/lapacke_ssbevd.c
  src/lapacke_ssbevd_work.c src/lapacke_ssbev_work.c src/lapacke_ssbevx.c
  src/lapacke_ssbevx_work.c src/lapacke_ssbgst.c src/lapacke_ssbgst_work.c
  src/lapacke_ssbgv.c src/lapacke_ssbgvd.c src/lapacke_ssbgvd_work.c
  src/lapacke_ssbgv_work.c src/lapacke_ssbgvx.c src/lapacke_ssbgvx_work.c
  src/lapacke_ssbtrd.c src/lapacke_ssbtrd_work.c src/lapacke_ssfrk.c
  src/lapacke_ssfrk_work.c src/lapacke_sspcon.c src/lapacke_sspcon_work.c
  src/lapacke_sspev.c src/lapacke_sspevd.c src/lapacke_sspevd_work.c
  src/lapacke_sspev_work.c src/lapacke_sspevx.c src/lapacke_sspevx_work.c
  src/lapacke_sspgst.c src/lapacke_sspgst_work.c src/lapacke_sspgv.c
  src/lapacke_sspgvd.c src/lapacke_sspgvd_work.c src/lapacke_sspgv_work.c
  src/lapacke_sspgvx.c src/lapacke_sspgvx_work.c src/lapacke_ssprfs.c
  src/lapacke_ssprfs_work.c src/lapacke_sspsv.c src/lapacke_sspsv_work.c
  src/lapacke_sspsvx.c src/lapacke_sspsvx_work.c src/lapacke_ssptrd.c
  src/lapacke_ssptrd_work.c src/lapacke_ssptrf.c src/lapacke_ssptrf_work.c
  src/lapacke_ssptri.c src/lapacke_ssptri_work.c src/lapacke_ssptrs.c
  src/lapacke_ssptrs_work.c src/lapacke_sstebz.c src/lapacke_sstebz_work.c
  src/lapacke_sstedc.c src/lapacke_sstedc_work.c src/lapacke_sstegr.c
  src/lapacke_sstegr_work.c src/lapacke_sstein.c src/lapacke_sstein_work.c
  src/lapacke_sstemr.c src/lapacke_sstemr_work.c src/lapacke_ssteqr.c
  src/lapacke_ssteqr_work.c src/lapacke_ssterf.c src/lapacke_ssterf_work.c
  src/lapacke_sstev.c src/lapacke_sstevd.c src/lapacke_sstevd_work.c
  src/lapacke_sstevr.c src/lapacke_sstevr_work.c src/lapacke_sstev_work.c
  src/lapacke_sstevx.c src/lapacke_sstevx_work.c src/lapacke_ssycon.c
  src/lapacke_ssycon_work.c src/lapacke_ssyequb.c src/lapacke_ssyequb_work.c
  src/lapacke_ssyev.c src/lapacke_ssyevd.c src/lapacke_ssyevd_work.c
  src/lapacke_ssyevr.c src/lapacke_ssyevr_work.c src/lapacke_ssyev_work.c
  src/lapacke_ssyevx.c src/lapacke_ssyevx_work.c src/lapacke_ssygst.c
  src/lapacke_ssygst_work.c src/lapacke_ssygv.c src/lapacke_ssygvd.c
  src/lapacke_ssygvd_work.c src/lapacke_ssygv_work.c src/lapacke_ssygvx.c
  src/lapacke_ssygvx_work.c src/lapacke_ssyrfs.c src/lapacke_ssyrfs_work.c
  src/lapacke_ssyrfsx.c src/lapacke_ssyrfsx_work.c src/lapacke_ssysv.c
  src/lapacke_ssysv_work.c src/lapacke_ssysvx.c src/lapacke_ssysvx_work.c
  src/lapacke_ssysvxx.c src/lapacke_ssysvxx_work.c src/lapacke_ssytrd.c
  src/lapacke_ssytrd_work.c src/lapacke_ssytrf.c src/lapacke_ssytrf_work.c
  src/lapacke_ssytri.c src/lapacke_ssytri_work.c src/lapacke_ssytrs.c
  src/lapacke_ssytrs_work.c src/lapacke_stbcon.c src/lapacke_stbcon_work.c
  src/lapacke_stbrfs.c src/lapacke_stbrfs_work.c src/lapacke_stbtrs.c
  src/lapacke_stbtrs_work.c src/lapacke_stfsm.c src/lapacke_stfsm_work.c
  src/lapacke_stftri.c src/lapacke_stftri_work.c src/lapacke_stfttp.c
  src/lapacke_stfttp_work.c src/lapacke_stfttr.c src/lapacke_stfttr_work.c
  src/lapacke_stgevc.c src/lapacke_stgevc_work.c src/lapacke_stgexc.c
  src/lapacke_stgexc_work.c src/lapacke_stgsen.c src/lapacke_stgsen_work.c
  src/lapacke_stgsja.c src/lapacke_stgsja_work.c src/lapacke_stgsna.c
  src/lapacke_stgsna_work.c src/lapacke_stgsyl.c src/lapacke_stgsyl_work.c
  src/lapacke_stpcon.c src/lapacke_stpcon_work.c src/lapacke_stprfs.c
  src/lapacke_stprfs_work.c src/lapacke_stptri.c src/lapacke_stptri_work.c
  src/lapacke_stptrs.c src/lapacke_stptrs_work.c src/lapacke_stpttf.c
  src/lapacke_stpttf_work.c src/lapacke_stpttr.c src/lapacke_stpttr_work.c
  src/lapacke_strcon.c src/lapacke_strcon_work.c src/lapacke_strevc.c
  src/lapacke_strevc_work.c src/lapacke_strexc.c src/lapacke_strexc_work.c
  src/lapacke_strrfs.c src/lapacke_strrfs_work.c src/lapacke_strsen.c
  src/lapacke_strsen_work.c src/lapacke_strsna.c src/lapacke_strsna_work.c
  src/lapacke_strsyl.c src/lapacke_strsyl_work.c src/lapacke_strtri.c
  src/lapacke_strtri_work.c src/lapacke_strtrs.c src/lapacke_strtrs_work.c
  src/lapacke_strttf.c src/lapacke_strttf_work.c src/lapacke_strttp.c
  src/lapacke_strttp_work.c src/lapacke_stzrzf.c src/lapacke_stzrzf_work.c
  src/lapacke_zbdsqr.c src/lapacke_zbdsqr_work.c src/lapacke_zcgesv.c
  src/lapacke_zcgesv_work.c src/lapacke_zcposv.c src/lapacke_zcposv_work.c
  src/lapacke_zgbbrd.c src/lapacke_zgbbrd_work.c src/lapacke_zgbcon.c
  src/lapacke_zgbcon_work.c src/lapacke_zgbequb.c src/lapacke_zgbequb_work.c
  src/lapacke_zgbequ.c src/lapacke_zgbequ_work.c src/lapacke_zgbrfs.c
  src/lapacke_zgbrfs_work.c src/lapacke_zgbrfsx.c src/lapacke_zgbrfsx_work.c
  src/lapacke_zgbsv.c src/lapacke_zgbsv_work.c src/lapacke_zgbsvx.c
  src/lapacke_zgbsvx_work.c src/lapacke_zgbsvxx.c src/lapacke_zgbsvxx_work.c
  src/lapacke_zgbtrf.c src/lapacke_zgbtrf_work.c src/lapacke_zgbtrs.c
  src/lapacke_zgbtrs_work.c src/lapacke_zgebak.c src/lapacke_zgebak_work.c
  src/lapacke_zgebal.c src/lapacke_zgebal_work.c src/lapacke_zgebrd.c
  src/lapacke_zgebrd_work.c src/lapacke_zgecon.c src/lapacke_zgecon_work.c
  src/lapacke_zgeequb.c src/lapacke_zgeequb_work.c src/lapacke_zgeequ.c
  src/lapacke_zgeequ_work.c src/lapacke_zgees.c src/lapacke_zgees_work.c
  src/lapacke_zgeesx.c src/lapacke_zgeesx_work.c src/lapacke_zgeev.c
  src/lapacke_zgeev_work.c src/lapacke_zgeevx.c src/lapacke_zgeevx_work.c
  src/lapacke_zgehrd.c src/lapacke_zgehrd_work.c src/lapacke_zgelqf.c
  src/lapacke_zgelqf_work.c src/lapacke_zgels.c src/lapacke_zgelsd.c
  src/lapacke_zgelsd_work.c src/lapacke_zgelss.c src/lapacke_zgelss_work.c
  src/lapacke_zgels_work.c src/lapacke_zgelsy.c src/lapacke_zgelsy_work.c
  src/lapacke_zgeqlf.c src/lapacke_zgeqlf_work.c src/lapacke_zgeqp3.c
  src/lapacke_zgeqp3_work.c src/lapacke_zgeqpf.c src/lapacke_zgeqpf_work.c
  src/lapacke_zgeqrf.c src/lapacke_zgeqrf_work.c src/lapacke_zgerfs.c
  src/lapacke_zgerfs_work.c src/lapacke_zgerfsx.c src/lapacke_zgerfsx_work.c
  src/lapacke_zgerqf.c src/lapacke_zgerqf_work.c src/lapacke_zgesdd.c
  src/lapacke_zgesdd_work.c src/lapacke_zgesv.c src/lapacke_zgesvd.c
  src/lapacke_zgesvd_work.c src/lapacke_zgesv_work.c src/lapacke_zgesvx.c
  src/lapacke_zgesvx_work.c src/lapacke_zgesvxx.c src/lapacke_zgesvxx_work.c
  src/lapacke_zgetrf.c src/lapacke_zgetrf_work.c src/lapacke_zgetri.c
  src/lapacke_zgetri_work.c src/lapacke_zgetrs.c src/lapacke_zgetrs_work.c
  src/lapacke_zggbak.c src/lapacke_zggbak_work.c src/lapacke_zggbal.c
  src/lapacke_zggbal_work.c src/lapacke_zgges.c src/lapacke_zgges_work.c
  src/lapacke_zggesx.c src/lapacke_zggesx_work.c src/lapacke_zggev.c
  src/lapacke_zggev_work.c src/lapacke_zggevx.c src/lapacke_zggevx_work.c
  src/lapacke_zggglm.c src/lapacke_zggglm_work.c src/lapacke_zgghrd.c
  src/lapacke_zgghrd_work.c src/lapacke_zgglse.c src/lapacke_zgglse_work.c
  src/lapacke_zggqrf.c src/lapacke_zggqrf_work.c src/lapacke_zggrqf.c
  src/lapacke_zggrqf_work.c src/lapacke_zggsvd.c src/lapacke_zggsvd_work.c
  src/lapacke_zggsvp.c src/lapacke_zggsvp_work.c src/lapacke_zgtcon.c
  src/lapacke_zgtcon_work.c src/lapacke_zgtrfs.c src/lapacke_zgtrfs_work.c
  src/lapacke_zgtsv.c src/lapacke_zgtsv_work.c src/lapacke_zgtsvx.c
  src/lapacke_zgtsvx_work.c src/lapacke_zgttrf.c src/lapacke_zgttrf_work.c
  src/lapacke_zgttrs.c src/lapacke_zgttrs_work.c src/lapacke_zhbev.c
  src/lapacke_zhbevd.c src/lapacke_zhbevd_work.c src/lapacke_zhbev_work.c
  src/lapacke_zhbevx.c src/lapacke_zhbevx_work.c src/lapacke_zhbgst.c
  src/lapacke_zhbgst_work.c src/lapacke_zhbgv.c src/lapacke_zhbgvd.c
  src/lapacke_zhbgvd_work.c src/lapacke_zhbgv_work.c src/lapacke_zhbgvx.c
  src/lapacke_zhbgvx_work.c src/lapacke_zhbtrd.c src/lapacke_zhbtrd_work.c
  src/lapacke_zhecon.c src/lapacke_zhecon_work.c src/lapacke_zheequb.c
  src/lapacke_zheequb_work.c src/lapacke_zheev.c src/lapacke_zheevd.c
  src/lapacke_zheevd_work.c src/lapacke_zheevr.c src/lapacke_zheevr_work.c
  src/lapacke_zheev_work.c src/lapacke_zheevx.c src/lapacke_zheevx_work.c
  src/lapacke_zhegst.c src/lapacke_zhegst_work.c src/lapacke_zhegv.c
  src/lapacke_zhegvd.c src/lapacke_zhegvd_work.c src/lapacke_zhegv_work.c
  src/lapacke_zhegvx.c src/lapacke_zhegvx_work.c src/lapacke_zherfs.c
  src/lapacke_zherfs_work.c src/lapacke_zherfsx.c src/lapacke_zherfsx_work.c
  src/lapacke_zhesv.c src/lapacke_zhesv_work.c src/lapacke_zhesvx.c
  src/lapacke_zhesvx_work.c src/lapacke_zhesvxx.c src/lapacke_zhesvxx_work.c
  src/lapacke_zhetrd.c src/lapacke_zhetrd_work.c src/lapacke_zhetrf.c
  src/lapacke_zhetrf_work.c src/lapacke_zhetri.c src/lapacke_zhetri_work.c
  src/lapacke_zhetrs.c src/lapacke_zhetrs_work.c src/lapacke_zhfrk.c
  src/lapacke_zhfrk_work.c src/lapacke_zhgeqz.c src/lapacke_zhgeqz_work.c
  src/lapacke_zhpcon.c src/lapacke_zhpcon_work.c src/lapacke_zhpev.c
  src/lapacke_zhpevd.c src/lapacke_zhpevd_work.c src/lapacke_zhpev_work.c
  src/lapacke_zhpevx.c src/lapacke_zhpevx_work.c src/lapacke_zhpgst.c
  src/lapacke_zhpgst_work.c src/lapacke_zhpgv.c src/lapacke_zhpgvd.c
  src/lapacke_zhpgvd_work.c src/lapacke_zhpgv_work.c src/lapacke_zhpgvx.c
  src/lapacke_zhpgvx_work.c src/lapacke_zhprfs.c src/lapacke_zhprfs_work.c
  src/lapacke_zhpsv.c src/lapacke_zhpsv_work.c src/lapacke_zhpsvx.c
  src/lapacke_zhpsvx_work.c src/lapacke_zhptrd.c src/lapacke_zhptrd_work.c
  src/lapacke_zhptrf.c src/lapacke_zhptrf_work.c src/lapacke_zhptri.c
  src/lapacke_zhptri_work.c src/lapacke_zhptrs.c src/lapacke_zhptrs_work.c
  src/lapacke_zhsein.c src/lapacke_zhsein_work.c src/lapacke_zhseqr.c
  src/lapacke_zhseqr_work.c src/lapacke_zpbcon.c src/lapacke_zpbcon_work.c
  src/lapacke_zpbequ.c src/lapacke_zpbequ_work.c src/lapacke_zpbrfs.c
  src/lapacke_zpbrfs_work.c src/lapacke_zpbstf.c src/lapacke_zpbstf_work.c
  src/lapacke_zpbsv.c src/lapacke_zpbsv_work.c src/lapacke_zpbsvx.c
  src/lapacke_zpbsvx_work.c src/lapacke_zpbtrf.c src/lapacke_zpbtrf_work.c
  src/lapacke_zpbtrs.c src/lapacke_zpbtrs_work.c src/lapacke_zpftrf.c
  src/lapacke_zpftrf_work.c src/lapacke_zpftri.c src/lapacke_zpftri_work.c
  src/lapacke_zpftrs.c src/lapacke_zpftrs_work.c src/lapacke_zpocon.c
  src/lapacke_zpocon_work.c src/lapacke_zpoequb.c src/lapacke_zpoequb_work.c
  src/lapacke_zpoequ.c src/lapacke_zpoequ_work.c src/lapacke_zporfs.c
  src/lapacke_zporfs_work.c src/lapacke_zporfsx.c src/lapacke_zporfsx_work.c
  src/lapacke_zposv.c src/lapacke_zposv_work.c src/lapacke_zposvx.c
  src/lapacke_zposvx_work.c src/lapacke_zposvxx.c src/lapacke_zposvxx_work.c
  src/lapacke_zpotrf.c src/lapacke_zpotrf_work.c src/lapacke_zpotri.c
  src/lapacke_zpotri_work.c src/lapacke_zpotrs.c src/lapacke_zpotrs_work.c
  src/lapacke_zppcon.c src/lapacke_zppcon_work.c src/lapacke_zppequ.c
  src/lapacke_zppequ_work.c src/lapacke_zpprfs.c src/lapacke_zpprfs_work.c
  src/lapacke_zppsv.c src/lapacke_zppsv_work.c src/lapacke_zppsvx.c
  src/lapacke_zppsvx_work.c src/lapacke_zpptrf.c src/lapacke_zpptrf_work.c
  src/lapacke_zpptri.c src/lapacke_zpptri_work.c src/lapacke_zpptrs.c
  src/lapacke_zpptrs_work.c src/lapacke_zpstrf.c src/lapacke_zpstrf_work.c
  src/lapacke_zptcon.c src/lapacke_zptcon_work.c src/lapacke_zpteqr.c
  src/lapacke_zpteqr_work.c src/lapacke_zptrfs.c src/lapacke_zptrfs_work.c
  src/lapacke_zptsv.c src/lapacke_zptsv_work.c src/lapacke_zptsvx.c
  src/lapacke_zptsvx_work.c src/lapacke_zpttrf.c src/lapacke_zpttrf_work.c
  src/lapacke_zpttrs.c src/lapacke_zpttrs_work.c src/lapacke_zspcon.c
  src/lapacke_zspcon_work.c src/lapacke_zsprfs.c src/lapacke_zsprfs_work.c
  src/lapacke_zspsv.c src/lapacke_zspsv_work.c src/lapacke_zspsvx.c
  src/lapacke_zspsvx_work.c src/lapacke_zsptrf.c src/lapacke_zsptrf_work.c
  src/lapacke_zsptri.c src/lapacke_zsptri_work.c src/lapacke_zsptrs.c
  src/lapacke_zsptrs_work.c src/lapacke_zstedc.c src/lapacke_zstedc_work.c
  src/lapacke_zstegr.c src/lapacke_zstegr_work.c src/lapacke_zstein.c
  src/lapacke_zstein_work.c src/lapacke_zstemr.c src/lapacke_zstemr_work.c
  src/lapacke_zsteqr.c src/lapacke_zsteqr_work.c src/lapacke_zsycon.c
  src/lapacke_zsycon_work.c src/lapacke_zsyequb.c src/lapacke_zsyequb_work.c
  src/lapacke_zsyrfs.c src/lapacke_zsyrfs_work.c src/lapacke_zsyrfsx.c
  src/lapacke_zsyrfsx_work.c src/lapacke_zsysv.c src/lapacke_zsysv_work.c
  src/lapacke_zsysvx.c src/lapacke_zsysvx_work.c src/lapacke_zsysvxx.c
  src/lapacke_zsysvxx_work.c src/lapacke_zsytrf.c src/lapacke_zsytrf_work.c
  src/lapacke_zsytri.c src/lapacke_zsytri_work.c src/lapacke_zsytrs.c
  src/lapacke_zsytrs_work.c src/lapacke_ztbcon.c src/lapacke_ztbcon_work.c
  src/lapacke_ztbrfs.c src/lapacke_ztbrfs_work.c src/lapacke_ztbtrs.c
  src/lapacke_ztbtrs_work.c src/lapacke_ztfsm.c src/lapacke_ztfsm_work.c
  src/lapacke_ztftri.c src/lapacke_ztftri_work.c src/lapacke_ztfttp.c
  src/lapacke_ztfttp_work.c src/lapacke_ztfttr.c src/lapacke_ztfttr_work.c
  src/lapacke_ztgevc.c src/lapacke_ztgevc_work.c src/lapacke_ztgexc.c
  src/lapacke_ztgexc_work.c src/lapacke_ztgsen.c src/lapacke_ztgsen_work.c
  src/lapacke_ztgsja.c src/lapacke_ztgsja_work.c src/lapacke_ztgsna.c
  src/lapacke_ztgsna_work.c src/lapacke_ztgsyl.c src/lapacke_ztgsyl_work.c
  src/lapacke_ztpcon.c src/lapacke_ztpcon_work.c src/lapacke_ztprfs.c
  src/lapacke_ztprfs_work.c src/lapacke_ztptri.c src/lapacke_ztptri_work.c
  src/lapacke_ztptrs.c src/lapacke_ztptrs_work.c src/lapacke_ztpttf.c
  src/lapacke_ztpttf_work.c src/lapacke_ztpttr.c src/lapacke_ztpttr_work.c
  src/lapacke_ztrcon.c src/lapacke_ztrcon_work.c src/lapacke_ztrevc.c
  src/lapacke_ztrevc_work.c src/lapacke_ztrexc.c src/lapacke_ztrexc_work.c
  src/lapacke_ztrrfs.c src/lapacke_ztrrfs_work.c src/lapacke_ztrsen.c
  src/lapacke_ztrsen_work.c src/lapacke_ztrsna.c src/lapacke_ztrsna_work.c
  src/lapacke_ztrsyl.c src/lapacke_ztrsyl_work.c src/lapacke_ztrtri.c
  src/lapacke_ztrtri_work.c src/lapacke_ztrtrs.c src/lapacke_ztrtrs_work.c
  src/lapacke_ztrttf.c src/lapacke_ztrttf_work.c src/lapacke_ztrttp.c
  src/lapacke_ztrttp_work.c src/lapacke_ztzrzf.c src/lapacke_ztzrzf_work.c
  src/lapacke_zungbr.c src/lapacke_zungbr_work.c src/lapacke_zunghr.c
  src/lapacke_zunghr_work.c src/lapacke_zunglq.c src/lapacke_zunglq_work.c
  src/lapacke_zungql.c src/lapacke_zungql_work.c src/lapacke_zungqr.c
  src/lapacke_zungqr_work.c src/lapacke_zungrq.c src/lapacke_zungrq_work.c
  src/lapacke_zungtr.c src/lapacke_zungtr_work.c src/lapacke_zunmbr.c
  src/lapacke_zunmbr_work.c src/lapacke_zunmhr.c src/lapacke_zunmhr_work.c
  src/lapacke_zunmlq.c src/lapacke_zunmlq_work.c src/lapacke_zunmql.c
  src/lapacke_zunmql_work.c src/lapacke_zunmqr.c src/lapacke_zunmqr_work.c
  src/lapacke_zunmrq.c src/lapacke_zunmrq_work.c src/lapacke_zunmrz.c
  src/lapacke_zunmrz_work.c src/lapacke_zunmtr.c src/lapacke_zunmtr_work.c
  src/lapacke_zupgtr.c src/lapacke_zupgtr_work.c src/lapacke_zupmtr.c
  src/lapacke_zupmtr_work.c)

set(LIBSWIFTNAV_REQUIRED 
  src/lapacke_dgelss.c src/lapacke_dgelss_work.c
  src/lapacke_dgelsy.c src/lapacke_dgelsy_work.c
  utils/lapacke_xerbla.c utils/lapacke_dge_trans.c utils/lapacke_d_nancheck.c utils/lapacke_dge_nancheck.c
)

#add_library(lapacke ${OBJ} ${OBJ_UTILS})
add_library(lapacke ${LIBSWIFTNAV_REQUIRED})
target_link_libraries(lapacke blas)
target_link_libraries(lapacke lapack)

