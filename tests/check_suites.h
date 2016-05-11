#ifndef CHECK_SUITES_H
#define CHECK_SUITES_H

Suite* almanac_suite(void);
Suite* dgnss_management_test_suite(void);
Suite* baseline_test_suite(void);
Suite* amb_kf_test_suite(void);
Suite* observation_test_suite(void);
Suite* pvt_test_suite(void);
Suite* coord_system_suite(void);
Suite* rtcm3_suite(void);
Suite* bits_suite(void);
Suite* memory_pool_suite(void);
Suite* edc_suite(void);
Suite* linear_algebra_suite(void);
Suite* sats_management_test_suite(void);
Suite* ambiguity_test_suite(void);
Suite* filter_utils_suite(void);
Suite* ephemeris_suite(void);
Suite* set_suite(void);
Suite* viterbi_suite(void);
Suite* time_test_suite(void);
Suite* ionosphere_suite(void);
Suite* signal_test_suite(void);
Suite* track_test_suite(void);
Suite* l2c_capability_test_suite(void);
Suite* cnav_test_suite(void);
Suite* glo_decoder_test_suite(void);
Suite* troposphere_suite(void);

#endif /* CHECK_SUITES_H */
