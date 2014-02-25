include_directories(${CMAKE_CURRENT_SOURCE_DIR}/INCLUDE)
# Oh dear, f2c what are you doing?
set(F2C_FAIL_FLAGS "-Wno-parentheses -Wno-unused-parameter -Wno-sign-compare -Wno-implicit-function-declaration -Wno-strict-prototypes -Wno-unused-but-set-variable -Wno-unused-variable -Wno-format")
add_subdirectory(F2CLIBS)
add_subdirectory(BLAS)
add_subdirectory(SRC)
#if (NOT CMAKE_CROSSCOMPILING)
  #add_subdirectory(TESTING)
#endif ()
