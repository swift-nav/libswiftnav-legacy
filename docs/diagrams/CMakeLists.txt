set(BUILD_DIAGRAMS ON)

find_package(LATEX)
if (NOT PDFLATEX_COMPILER)
  message(WARNING
    "Couldn't find pdfLaTeX, diagrams will not be included in the "
    "documentation."
  )
  set(BUILD_DIAGRAMS OFF)
endif (NOT PDFLATEX_COMPILER)

find_program(IMAGEMAGICK_CONVERT convert
  DOC "The convert program (part of ImageMagick)."
)
if (NOT IMAGEMAGICK_CONVERT)
  message(WARNING
    "Couldn't find the convert program (part of ImageMagick), "
    "diagrams will not be included in the documentation."
  )
  set(BUILD_DIAGRAMS OFF)
endif (NOT IMAGEMAGICK_CONVERT)

function(build_diagrams out_var)
  set(result)
  foreach(in_f ${ARGN})

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${in_f}.png
      COMMENT "Building LaTeX diagram ${in_f}.tex"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND ${PDFLATEX_COMPILER}
        #-interaction=batchmode
        -output-directory ${CMAKE_CURRENT_BINARY_DIR}
        ${in_f}.tex
      COMMAND ${IMAGEMAGICK_CONVERT}
        -density 120 -trim -bordercolor White -border 5
        ${CMAKE_CURRENT_BINARY_DIR}/${in_f}.pdf
        ${CMAKE_CURRENT_BINARY_DIR}/${in_f}.png
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${in_f}.tex
    )

    list(APPEND result ${CMAKE_CURRENT_BINARY_DIR}/${in_f}.png)
  endforeach()
  set(${out_var} "${result}" PARENT_SCOPE)
endfunction()

if (BUILD_DIAGRAMS)
  build_diagrams(DIAGRAM_PNGS
    costas_loop
    2nd_order_dpll
    1st_order_loop_filter
    wgsecef2azel
  )

  add_custom_target(diagrams DEPENDS ${DIAGRAM_PNGS})

  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
    ${DIAGRAM_PNGS}
  )
else (BUILD_DIAGRAMS)
  add_custom_target(diagrams)
endif (BUILD_DIAGRAMS)

