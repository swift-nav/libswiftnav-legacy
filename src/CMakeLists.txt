if (NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif (NOT DEFINED BUILD_SHARED_LIBS)

file(GLOB libswiftnav_HEADERS "${PROJECT_SOURCE_DIR}/include/*.h")

include_directories("${PROJECT_SOURCE_DIR}/include")

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
)

if (CMAKE_CROSSCOMPILING)
  set(LIB_SUFFIX "-${CMAKE_SYSTEM_PROCESSOR}")
endif (CMAKE_CROSSCOMPILING)

add_library(swiftnav-static${LIB_SUFFIX} STATIC ${libswiftnav_SRCS})
install(TARGETS swiftnav-static${LIB_SUFFIX} DESTINATION lib)

if(BUILD_SHARED_LIBS)
  add_library(swiftnav${LIB_SUFFIX} SHARED ${libswiftnav_SRCS})
  install(TARGETS swiftnav${LIB_SUFFIX} DESTINATION lib)
else(BUILD_SHARED_LIBS)
  message(STATUS "Not building shared libraries")
endif(BUILD_SHARED_LIBS)

install(FILES ${libswiftnav_HEADERS} DESTINATION include/libswiftnav)

