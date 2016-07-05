##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-annotate/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
set ( SRC_FILES  main.cpp )
                                           
add_executable ( cfm-annotate ${SRC_FILES} )
target_link_libraries ( cfm-annotate cfm-code ${REQUIRED_LIBS} )

if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	install ( TARGETS cfm-annotate 
          	  DESTINATION bin )
else()
	install ( TARGETS cfm-annotate 
          	  DESTINATION ${CFM_OUTPUT_DIR} )
endif()
          
if (WIN32)
    install ( FILES ${LPSOLVE_INCLUDE_DIR}/lpsolve55.dll
              DESTINATION ${CFM_OUTPUT_DIR} )
endif (WIN32)              
          
