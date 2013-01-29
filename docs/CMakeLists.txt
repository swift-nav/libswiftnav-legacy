#option(BUILD_DOCUMENTATION "Use Doxygen to create API documentation" ON)
set(BUILD_DOCUMENTATION ON)

find_package(Doxygen)
if (NOT DOXYGEN_FOUND)
  set(BUILD_DOCUMENTATION OFF)
  message(WARNING
    "Doxygen not found, the documentation will not be built."
  )
endif (NOT DOXYGEN_FOUND)

if (BUILD_DOCUMENTATION)

  add_subdirectory(diagrams)

  configure_file(Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    @ONLY
  )

  add_custom_target(docs
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMENT "Building Doxygen documentation"
    DEPENDS diagrams Doxyfile.in
  )

  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
		${CMAKE_CURRENT_BINARY_DIR}/html
  )

endif (BUILD_DOCUMENTATION)

