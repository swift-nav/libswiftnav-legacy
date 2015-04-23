if (NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif (NOT DEFINED BUILD_SHARED_LIBS)

set(CMAKE_C_FLAGS "-Wmissing-prototypes ${CMAKE_C_FLAGS}")

file(GLOB libswiftnav_HEADERS "${PROJECT_SOURCE_DIR}/include/libswiftnav/*.h")

include_directories("${PROJECT_SOURCE_DIR}/CBLAS/include")
include_directories("${PROJECT_SOURCE_DIR}/clapack-3.2.1-CMAKE/INCLUDE")

include_directories("${PROJECT_SOURCE_DIR}/include/libswiftnav")

set_source_files_properties(${plover_SRCS} PROPERTIES GENERATED TRUE)
set_source_files_properties(${plover_HDRS} PROPERTIES GENERATED TRUE)

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
  lambda.c
  amb_kf.c
  baseline.c
  single_diff.c
  memory_pool.c
  dgnss_management.c
  sats_management.c
  ambiguity_test.c
  printing_utils.c
  ${plover_SRCS}
)

add_library(swiftnav-static STATIC ${libswiftnav_SRCS})
add_dependencies(swiftnav-static generate)
target_link_libraries(swiftnav-static cblas)
target_link_libraries(swiftnav-static lapack)
install(TARGETS swiftnav-static DESTINATION lib${LIB_SUFFIX})

if(BUILD_SHARED_LIBS)
  add_library(swiftnav SHARED ${libswiftnav_SRCS})
  add_dependencies(swiftnav generate)
  target_link_libraries(swiftnav cblas)
  target_link_libraries(swiftnav lapack)
  install(TARGETS swiftnav DESTINATION lib${LIB_SUFFIX})
else(BUILD_SHARED_LIBS)
  message(STATUS "Not building shared libraries")
endif(BUILD_SHARED_LIBS)

install(FILES ${libswiftnav_HEADERS} DESTINATION include/libswiftnav)

