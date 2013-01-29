find_package(Check REQUIRED)
include_directories(${CHECK_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CHECK_LIBRARIES} swiftnav)

include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(test_libswiftnav
  check_main.c
  check_utils.c

  check_coord_system.c
  check_linear_algebra.c
)

target_link_libraries(test_libswiftnav ${LIBS})
add_custom_command(
  TARGET test_libswiftnav POST_BUILD
  COMMENT "Running unit tests"
  COMMAND test_libswiftnav
)

