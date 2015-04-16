# Working dir
set (CMAKE_PLOVER_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")

# Find default cabal
find_program(CMAKE_cabal NAMES cabal)

if (DEFINED ENV{CABAL})
  message("Saw cabal flag: $ENV{CABAL}")
  set (CMAKE_cabal $ENV{CABAL})
else ()
  message("Environment variable 'CABAL' not set. Defaulting to ${CMAKE_cabal}.")
endif ()

if (CMAKE_cabal)
  message("\nFound Haskell tools. Adding code generation step (to be run if Haskell subproject is modified).\n")
  add_custom_command(
    OUTPUT "${CMAKE_PLOVER_DIR}/generated.h" "${CMAKE_PLOVER_DIR}/generated.c"
    COMMAND ${CMAKE_cabal} run "${CMAKE_PLOVER_DIR}/generated"
    WORKING_DIRECTORY ${CMAKE_PLOVER_DIR}
    DEPENDS Main.hs)
else ()
  message("\nNo Haskell toolchain (cabal) found. Code generation with Plover is disabled.\n")
endif ()

# Main project depends on generated code
add_custom_target(generate DEPENDS "${CMAKE_PLOVER_DIR}/generated.c")
