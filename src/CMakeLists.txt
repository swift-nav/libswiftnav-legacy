if (NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif (NOT DEFINED BUILD_SHARED_LIBS)

set(CMAKE_C_FLAGS "-Wmissing-prototypes ${CMAKE_C_FLAGS}")

file(GLOB libswiftnav_HEADERS "${PROJECT_SOURCE_DIR}/include/libswiftnav/*.h")

include_directories("${PROJECT_SOURCE_DIR}/CBLAS/include")
include_directories("${PROJECT_SOURCE_DIR}/clapack-3.2.1-CMAKE/INCLUDE")
include_directories("${PROJECT_SOURCE_DIR}/libfec/include")

include_directories("${PROJECT_SOURCE_DIR}/include")

set_source_files_properties(${plover_SRCS} PROPERTIES GENERATED TRUE)
set_source_files_properties(${plover_HDRS} PROPERTIES GENERATED TRUE)

set(libswiftnav_SRCS
  logging.c
  ephemeris.c
  nav_msg.c
  pvt.c
  troposphere.c
  track.c
  correlate.c
  coord_system.c
  linear_algebra.c
  prns.c
  almanac.c
  time.c
  edc.c
  rtcm3.c
  bits.c
  lambda.c
  amb_kf.c
  baseline.c
  observation.c
  set.c
  memory_pool.c
  dgnss_management.c
  sats_management.c
  ambiguity_test.c
  printing_utils.c
  filter_utils.c
  signal.c
  ionosphere.c
  bit_sync.c
  cnav_msg.c
  ${plover_SRCS}

  CACHE INTERNAL ""
)

add_library(swiftnav-static STATIC ${libswiftnav_SRCS})
add_dependencies(swiftnav-static generate)
target_link_libraries(swiftnav-static cblas)
target_link_libraries(swiftnav-static lapack)
target_link_libraries(swiftnav-static fec)
install(TARGETS swiftnav-static DESTINATION lib${LIB_SUFFIX})

if(BUILD_SHARED_LIBS)
  add_library(swiftnav SHARED ${libswiftnav_SRCS})
  add_dependencies(swiftnav generate)
  target_link_libraries(swiftnav cblas)
  target_link_libraries(swiftnav lapack)
  target_link_libraries(swiftnav fec)
  install(TARGETS swiftnav DESTINATION lib${LIB_SUFFIX})
else(BUILD_SHARED_LIBS)
  message(STATUS "Not building shared libraries")
endif(BUILD_SHARED_LIBS)

install(FILES ${libswiftnav_HEADERS} DESTINATION include/libswiftnav)

MESSAGE("${PROJECT_SOURCE_DIR}")

add_custom_command(
  OUTPUT checkpatch
  DEPENDS ${libswiftnav_SRCS}
  VERBATIM
  COMMAND "${PROJECT_SOURCE_DIR}/checks/checkoutput.sh")

add_custom_target(checkpatch DEPENDS checkpatch)
