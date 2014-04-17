# Disable code coverage for this subdirectory
set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")

# Build libraries with PIC flag as they will be linked into shared objects
if (NOT CMAKE_CROSSCOMPILING)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CBLAS_FAIL_FLAGS "-Wno-char-subscripts")
add_subdirectory(src)
