include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CBLAS_FAIL_FLAGS "-Wno-char-subscripts")
add_subdirectory(src)
