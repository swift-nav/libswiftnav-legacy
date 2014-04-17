# Disable code coverage for this subdirectory
set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}")
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")

# Build libraries with PIC flag as they will be linked into shared objects
if (NOT CMAKE_CROSSCOMPILING)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/INCLUDE)
# Oh dear, f2c what are you doing?
set(F2C_FAIL_FLAGS "-Wno-parentheses -Wno-unused-parameter -Wno-sign-compare -Wno-implicit-function-declaration -Wno-strict-prototypes -Wno-unused-but-set-variable -Wno-unused-variable -Wno-format -Wno-format-security -Wno-uninitialized -Wno-self-assign")
add_subdirectory(F2CLIBS)
add_subdirectory(BLAS)
add_subdirectory(SRC)
#if (NOT CMAKE_CROSSCOMPILING)
  #add_subdirectory(TESTING)
#endif ()
