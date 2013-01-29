project(libswiftnav)

# Set project version using Git tag and hash.
execute_process(
  COMMAND "${PROJECT_SOURCE_DIR}/scripts/gitversion"
  OUTPUT_VARIABLE VERSION
)

cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Some compiler options used globally
set(CMAKE_C_FLAGS "-Wall -Wextra -Werror -std=gnu99 ${CMAKE_C_FLAGS}")

add_subdirectory(src)
add_subdirectory(docs)
add_subdirectory(tests)

