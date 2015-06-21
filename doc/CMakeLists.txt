FIND_PROGRAM( DOXYGEN_CMD doxygen )

# Configure the Doxyfile
CONFIGURE_FILE( Doxyfile.in Doxyfile )

# Can type "make doc" to generate documentation.
ADD_CUSTOM_TARGET( doc )

IF( DOXYGEN_CMD )
   ADD_CUSTOM_TARGET(
      doc-doxygen
      COMMAND ${DOXYGEN_CMD} Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   )
   ADD_DEPENDENCIES(doc doc-doxygen)
   INSTALL( DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html" DESTINATION doc/nvl )
ELSE()
   MESSAGE( STATUS "You do not have Doxygen. Documentation will not be available." )
ENDIF()