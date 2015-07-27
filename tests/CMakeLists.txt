if (CMAKE_CROSSCOMPILING)
  message(STATUS "Skipping unit tests, cross compiling")
else (CMAKE_CROSSCOMPILING)

  find_package(Check)
  if (NOT CHECK_FOUND)
    message(STATUS "Skipping unit tests, Check library not found!")
  else (NOT CHECK_FOUND)

    include_directories("${PROJECT_SOURCE_DIR}/CBLAS/include")
    include_directories("${PROJECT_SOURCE_DIR}/clapack-3.2.1-CMAKE/INCLUDE")

    include_directories(${CHECK_INCLUDE_DIRS})
    set(TEST_LIBS ${TEST_LIBS} ${CHECK_LIBRARIES} pthread swiftnav lapack cblas m fec)
    # Check needs to be linked against Librt on Linux
    if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
      set(TEST_LIBS ${TEST_LIBS} rt)
    endif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

    include_directories("${PROJECT_SOURCE_DIR}/include/libswiftnav")
    include_directories("${PROJECT_SOURCE_DIR}/src")
    include_directories("${PROJECT_SOURCE_DIR}/libfec/include")

    add_executable(test_libswiftnav
      check_main.c
      check_utils.c

      check_dgnss_management.c
      check_baseline.c
      check_amb_kf.c
      check_observation.c
      check_pvt.c
      check_edc.c
      check_bits.c
      check_memory_pool.c
      check_rtcm3.c
      check_coord_system.c
      check_linear_algebra.c
      check_sats_management.c
      check_ambiguity_test.c
      check_filter_utils.c
      check_ephemeris.c
      check_set.c
      check_viterbi.c
      check_gpstime.c
    )

    target_link_libraries(test_libswiftnav ${TEST_LIBS})
    add_custom_command(
      TARGET test_libswiftnav POST_BUILD
      COMMENT "Running unit tests"
      COMMAND test_libswiftnav
    )
    file(COPY waas_check.bin DESTINATION ${CMAKE_BINARY_DIR}/tests)
    file(COPY waas_data.bin DESTINATION ${CMAKE_BINARY_DIR}/tests)

  endif (NOT CHECK_FOUND)
endif (CMAKE_CROSSCOMPILING)
