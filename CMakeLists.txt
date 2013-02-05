project(libswiftnav)

# Set project version using Git tag and hash.
execute_process(
  COMMAND git describe --dirty --tags --always
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  RESULT_VARIABLE GIT_VERSION_FOUND
  ERROR_QUIET
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (GIT_VERSION_FOUND)
  set(VERSION "unknown")
else (GIT_VERSION_FOUND)
  set(VERSION ${GIT_VERSION})
endif (GIT_VERSION_FOUND)

# Set project version explicitly for release tarballs.
#set(VERSION foo)

message(STATUS "libswiftnav version: ${VERSION}")

cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Some compiler options used globally
set(CMAKE_C_FLAGS "-Wall -Wextra -Werror -std=gnu99 ${CMAKE_C_FLAGS}")

if (NOT CMAKE_CROSSCOMPILING)
  # Detect and use optimised compiler flags for the host architecture,
  # this is specific to x86 family CPUs.
  include(cmake/OptimizeForArchitecture.cmake)
  OptimizeForArchitecture()
  message(STATUS "Adding architecture flags: ${Vc_ARCHITECTURE_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${Vc_ARCHITECTURE_FLAGS}")
endif ()

add_subdirectory(src)
add_subdirectory(docs)
add_subdirectory(tests)

