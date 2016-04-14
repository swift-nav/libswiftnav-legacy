# Working dir
set (PLOVER_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")

# Output dirs
set (PLOVER_GENERATED_C_DIR ${CMAKE_SOURCE_DIR}/src/plover CACHE INTERNAL "")
file(MAKE_DIRECTORY ${PLOVER_GENERATED_C_DIR})
set (PLOVER_GENERATED_H_DIR ${CMAKE_SOURCE_DIR}/include/libswiftnav/plover CACHE INTERNAL "")
file(MAKE_DIRECTORY ${PLOVER_GENERATED_H_DIR})

# Find default cabal
find_program(CMAKE_CABAL NAMES cabal)
# Find plover executable
find_program(CMAKE_PLOVER NAMES plover)

if (DEFINED ENV{CABAL})
  message("Saw CABAL flag: $ENV{CABAL}")
  set (CMAKE_CABAL $ENV{CABAL})
else ()
  message("Environment variable 'CABAL' not set. Defaulting to ${CMAKE_CABAL}.")
endif ()
if (DEFINED ENV{PLOVER})
  message("Saw PLOVER flag. Using $ENV{PLOVER}")
  set (CMAKE_PLOVER $ENV{PLOVER})
else ()
  message("Environment variable 'PLOVER' not set. Defaulting to ${CMAKE_PLOVER}.")
endif ()

# Haskell files with macro plover definitions
# TODO implement build step
set (PLOVER_HS_SRCS
)

# For each unit in PLOVER_UNITS:
#   Expects ${unit}.plv in PLOVER_DIR
#   Generates ${unit}.c in PLOVER_GENERATED_C_DIR
#         and ${unit}.h in PLOVER_GENERATED_H_DIR
set (PLOVER_UNITS
  prelude
  qr
)

foreach(unit ${PLOVER_UNITS})
  list(APPEND PLOVER_CFILES "${PLOVER_GENERATED_C_DIR}/${unit}.c")
  list(APPEND PLOVER_HFILES "${PLOVER_GENERATED_H_DIR}/${unit}.h")
  list(APPEND PLOVER_SRC    "${PLOVER_DIR}/${unit}.plv")
endforeach(unit)

# Make available to src/ CMakeLists
set(PLOVER_CFILES ${PLOVER_CFILES} PARENT_SCOPE)
set(PLOVER_HFILES ${PLOVER_HFILES} PARENT_SCOPE)

# Main project depends on generate. generate depends on generated files.
add_custom_target(generate DEPENDS ${PLOVER_CFILES} ${PLOVER_HFILES})

if (CMAKE_CABAL AND CMAKE_PLOVER)
  message("\nFound Plover tools.")
  # Don't add step if no UNITS are specified
  if (PLOVER_CFILES OR PLOVER_HFILES)
    message("Adding code generation steps, one for each plv unit.\n\
Will run if plover subproject is modified and after make clean or clean-generate.")
    foreach(unit ${PLOVER_UNITS})
      add_custom_command(
        OUTPUT "${PLOVER_GENERATED_C_DIR}/${unit}.c" "${PLOVER_GENERATED_H_DIR}/${unit}.h"
        COMMAND ${CMAKE_PLOVER}
          -I "${PLOVER_DIR}"
          --cdir "${PLOVER_GENERATED_C_DIR}" --hdir "${PLOVER_GENERATED_H_DIR}"
          --libprefix "libswiftnav/plover"
          "${PLOVER_DIR}/${unit}.plv"
        WORKING_DIRECTORY ${PLOVER_DIR}
        DEPENDS "${PLOVER_DIR}/${unit}.plv")
    endforeach(unit)
  endif ()
  message("")

else ()
  message("\nNo Plover toolchain found. Code generation is disabled.\n")
endif ()

add_custom_target(clean-generate
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${PLOVER_DIR}/dist
  COMMAND ${CMAKE_COMMAND} -E remove ${PLOVER_GENERATED_C_DIR}/*
  COMMAND ${CMAKE_COMMAND} -E remove ${PLOVER_GENERATED_H_DIR}/*
)
