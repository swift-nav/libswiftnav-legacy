add_custom_command (
  OUTPUT "${PROJECT_SOURCE_DIR}/plover/generated.c" "${PROJECT_SOURCE_DIR}/plover/generated.h"
  COMMAND cabal run
  DEPENDS "${PROJECT_SOURCE_DIR}/plover/Main.hs"
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/plover/"
)

add_custom_target(generate DEPENDS "${PROJECT_SOURCE_DIR}/plover/generated.c" "${PROJECT_SOURCE_DIR}/plover/generated.h")
add_dependencies(swiftnav-static generate)
if(BUILD_SHARED_LIBS)
  add_dependencies(swiftnav generate)
endif(BUILD_SHARED_LIBS)
