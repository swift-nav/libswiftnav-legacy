if (NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif (NOT DEFINED BUILD_SHARED_LIBS)

set(CMAKE_C_FLAGS "-Wmissing-prototypes ${CMAKE_C_FLAGS}")

file(GLOB libfec_HEADERS "${PROJECT_SOURCE_DIR}/libfec/include/*.h")

set(viterbi27.c)
add_library(fec viterbi27.c)

install(TARGETS fec DESTINATION lib${LIB_SUFFIX})

install(FILES ${libfec_HEADERS} DESTINATION include/libswiftnav)
