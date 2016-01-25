if (NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif (NOT DEFINED BUILD_SHARED_LIBS)

set(CMAKE_C_FLAGS "-Wmissing-prototypes ${CMAKE_C_FLAGS}")

file(GLOB libl2cbitstream_HEADERS
"${PROJECT_SOURCE_DIR}/tests/data/l2cbitstream/libl2cbitstream/l2cbitstream.h")

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories(".")

set(libl2cbitstream_SRCS
  l2cbitstream.c

  CACHE INTERNAL ""
)

add_library(l2cbitstream-static STATIC ${libl2cbitstream_SRCS})
add_dependencies(l2cbitstream-static generate)
target_link_libraries(l2cbitstream-static swiftnav)
install(TARGETS l2cbitstream-static DESTINATION lib${LIB_SUFFIX})

if(BUILD_SHARED_LIBS)
  add_library(l2cbitstream SHARED ${libl2cbitstream_SRCS})
  add_dependencies(l2cbitstream generate)
  target_link_libraries(l2cbitstream swiftnav)
  install(TARGETS l2cbitstream DESTINATION lib${LIB_SUFFIX})
else(BUILD_SHARED_LIBS)
  message(STATUS "Not building shared libraries")
endif(BUILD_SHARED_LIBS)

install(FILES ${libl2cbitstream_HEADERS} DESTINATION include/libl2cbitstream)

MESSAGE("${PROJECT_SOURCE_DIR}")

