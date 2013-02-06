if (CMAKE_CROSSCOMPILING)
  message(STATUS "Skipping unit tests, cross compiling")
else (CMAKE_CROSSCOMPILING)

  find_package(Check)
  if (NOT CHECK_FOUND)
    message(STATUS "Skipping unit tests, Check library not found!")
  else (NOT CHECK_FOUND)

    include_directories(${CHECK_INCLUDE_DIRS})
    set(TEST_LIBS ${TEST_LIBS} ${CHECK_LIBRARIES} swiftnav m)

    include_directories("${PROJECT_SOURCE_DIR}/include")

    add_executable(test_libswiftnav
      check_main.c
      check_utils.c

      check_coord_system.c
      check_linear_algebra.c
    )

    target_link_libraries(test_libswiftnav ${TEST_LIBS})
    add_custom_command(
      TARGET test_libswiftnav POST_BUILD
      COMMENT "Running unit tests"
      COMMAND test_libswiftnav
    )

  endif (NOT CHECK_FOUND)
endif (CMAKE_CROSSCOMPILING)

