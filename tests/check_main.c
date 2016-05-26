#include <stdlib.h>
#include <check.h>

#include "check_suites.h"

int main(void)
{
  int number_failed;

  Suite *s = edc_suite();

  SRunner *sr = srunner_create(s);
  srunner_set_xml(sr, "test_results.xml");

  srunner_add_suite(sr, almanac_suite());
  srunner_add_suite(sr, dgnss_management_test_suite());
  srunner_add_suite(sr, baseline_test_suite());
  srunner_add_suite(sr, amb_kf_test_suite());
  srunner_add_suite(sr, observation_test_suite());
  srunner_add_suite(sr, pvt_test_suite());
  srunner_add_suite(sr, sats_management_test_suite());
  srunner_add_suite(sr, ambiguity_test_suite());
  srunner_add_suite(sr, rtcm3_suite());
  srunner_add_suite(sr, bits_suite());
  srunner_add_suite(sr, memory_pool_suite());
  srunner_add_suite(sr, coord_system_suite());
  srunner_add_suite(sr, linear_algebra_suite());
  srunner_add_suite(sr, filter_utils_suite());
  srunner_add_suite(sr, ephemeris_suite());
  srunner_add_suite(sr, set_suite());
  srunner_add_suite(sr, viterbi_suite());
  srunner_add_suite(sr, time_test_suite());
  srunner_add_suite(sr, ionosphere_suite());
  srunner_add_suite(sr, signal_test_suite());
  srunner_add_suite(sr, track_test_suite());
  srunner_add_suite(sr, l2c_capability_test_suite());
  srunner_add_suite(sr, cnav_test_suite());
  srunner_add_suite(sr, glo_decoder_test_suite());
  srunner_add_suite(sr, troposphere_suite());
  srunner_add_suite(sr, correlator_suite());
  srunner_add_suite(sr, counter_checker_suite());
  srunner_add_suite(sr, prns_test_suite());

  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
