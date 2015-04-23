# Working dir
set (PLOVER_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")

# Output dirs
set (PLOVER_GENERATED_C_DIR ${CMAKE_SOURCE_DIR}/src/plover CACHE INTERNAL "")
file(MAKE_DIRECTORY ${PLOVER_GENERATED_C_DIR})
set (PLOVER_GENERATED_H_DIR ${CMAKE_SOURCE_DIR}/include/libswiftnav/plover CACHE INTERNAL "")
file(MAKE_DIRECTORY ${PLOVER_GENERATED_H_DIR})

# Find default cabal
find_program(CMAKE_cabal NAMES cabal)

if (DEFINED ENV{CABAL})
  message("Saw cabal flag: $ENV{CABAL}")
  set (CMAKE_cabal $ENV{CABAL})
else ()
  message("Environment variable 'CABAL' not set. Defaulting to ${CMAKE_cabal}.")
endif ()

set (plover_HS_SRCS
  Main.hs
  AmbiguityTest.hs
)

set (plover_SRCS
  ${PLOVER_GENERATED_C_DIR}/ambiguity_test.c
  CACHE INTERNAL ""
)

set (plover_HDRS
  ${PLOVER_GENERATED_H_DIR}/ambiguity_test.h
  CACHE INTERNAL ""
)

if (CMAKE_cabal)
  message("\nFound Haskell tools. Adding code generation step (to be run if Haskell subproject is modified).\n")
  add_custom_command(
    OUTPUT ${plover_SRCS} ${plover_HDRS}
    COMMAND ${CMAKE_cabal} run "${PLOVER_GENERATED_C_DIR}" "${PLOVER_GENERATED_H_DIR}"
    WORKING_DIRECTORY ${PLOVER_DIR}
    DEPENDS ${plover_HS_SRCS})
else ()
  message("\nNo Haskell toolchain (cabal) found. Code generation with Plover is disabled.\n")
endif ()

add_custom_target(clean_generate
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${PLOVER_DIR}/dist
  COMMAND ${CMAKE_COMMAND} -E remove ${PLOVER_GENERATED_C_DIR}/*
  COMMAND ${CMAKE_COMMAND} -E remove ${PLOVER_GENERATED_H_DIR}/*
)

# Main project depends on generated code
add_custom_target(generate DEPENDS ${plover_SRCS} ${plover_HDRS})

