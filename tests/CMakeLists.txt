if (CMAKE_CROSSCOMPILING)
  message(STATUS "Skipping unit tests, cross compiling")
else (CMAKE_CROSSCOMPILING)

  find_package(Check)
  if (NOT CHECK_FOUND)
    message(STATUS "Skipping unit tests, Check library not found!")
  else (NOT CHECK_FOUND)

    include_directories("${PROJECT_SOURCE_DIR}/CBLAS/include")
    include_directories("${PROJECT_SOURCE_DIR}/clapack-3.2.1-CMAKE/INCLUDE")
    include_directories("${PROJECT_SOURCE_DIR}/lapacke/include")

    include_directories(${CHECK_INCLUDE_DIRS})
    set(TEST_LIBS ${TEST_LIBS} ${CHECK_LIBRARIES} pthread swiftnav lapacke cblas m)

    include_directories("${PROJECT_SOURCE_DIR}/include/libswiftnav")

    add_executable(test_libswiftnav
      check_main.c
      check_utils.c

      check_float_kf.c
      check_edc.c
      check_bits.c
      check_memory_pool.c
      check_sbp.c
      check_rtcm3.c
      check_coord_system.c
      check_linear_algebra.c
      check_ambiguity_test.c
    )

    target_link_libraries(test_libswiftnav ${TEST_LIBS})
    add_custom_command(
      TARGET test_libswiftnav POST_BUILD
      COMMENT "Running unit tests"
      COMMAND test_libswiftnav
    )

  endif (NOT CHECK_FOUND)
endif (CMAKE_CROSSCOMPILING)

