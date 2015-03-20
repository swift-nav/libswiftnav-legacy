add_custom_command (
  OUTPUT "${PROJECT_SOURCE_DIR}/plover/plover.c" "${PROJECT_SOURCE_DIR}/plover/plover.h"
  COMMAND cabal run
  DEPENDS "${PROJECT_SOURCE_DIR}/plover/Main.hs"
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/plover/"
)

add_custom_target(PloverGeneration DEPENDS "${PROJECT_SOURCE_DIR}/plover/plover.c")
add_dependencies(swiftnav-static PloverGeneration)
add_dependencies(swiftnav PloverGeneration)
