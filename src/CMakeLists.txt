if (NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif (NOT DEFINED BUILD_SHARED_LIBS)

file(GLOB libswiftnav_HEADERS "${PROJECT_SOURCE_DIR}/include/libswiftnav/*.h")

include_directories("${PROJECT_SOURCE_DIR}/CBLAS/include")
include_directories("${PROJECT_SOURCE_DIR}/clapack-3.2.1-CMAKE/INCLUDE")
include_directories("${PROJECT_SOURCE_DIR}/lapacke/include")

include_directories("${PROJECT_SOURCE_DIR}/include/libswiftnav")

set(libswiftnav_SRCS
  ephemeris.c
  nav_msg.c
  pvt.c
  tropo.c
  track.c
  correlate.c
  coord_system.c
  linear_algebra.c
  prns.c
  almanac.c
  gpstime.c
  edc.c
  rtcm3.c
  bits.c
  sbp.c
  lambda.c
  float_kf.c
  amb_kf.c
  stupid_filter.c
  sbp_utils.c
  single_diff.c
  memory_pool.c
  dgnss_management.c
  sats_management.c
  ambiguity_test.c
)

add_library(swiftnav-static STATIC ${libswiftnav_SRCS})
target_link_libraries(swiftnav-static cblas)
target_link_libraries(swiftnav-static lapacke)
install(TARGETS swiftnav-static DESTINATION lib${LIB_SUFFIX})

if(BUILD_SHARED_LIBS)
  add_library(swiftnav SHARED ${libswiftnav_SRCS})
  target_link_libraries(swiftnav cblas)
  target_link_libraries(swiftnav lapacke)
  install(TARGETS swiftnav DESTINATION lib${LIB_SUFFIX})
else(BUILD_SHARED_LIBS)
  message(STATUS "Not building shared libraries")
endif(BUILD_SHARED_LIBS)

install(FILES ${libswiftnav_HEADERS} DESTINATION include/libswiftnav)

