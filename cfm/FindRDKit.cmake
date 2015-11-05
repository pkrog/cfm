# FindRDKit.cmake
# Placed in the public domain by NextMove Software in 2013
# Modified by Felicity Allen, August 2013
# Try to find RDKit headers and libraries
# Defines:
#
#  RDKIT_FOUND - system has RDKit
#  RDKIT_INCLUDE_DIR - the RDKit include directory
#  RDKIT_INCLUDE_EXT_DIR - the RDKit external directory when including Inchi support
#  RDKIT_LIBRARIES - Link these to use RDKit

if(RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
  # in cache already or user-specified
  set(RDKIT_FOUND TRUE)

else()

  if(NOT RDKIT_INCLUDE_DIR)
    if(WIN32)
      find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
        PATHS
        ${RDKIT_DIR}\\Code
        $ENV{RDKIT_INCLUDE_DIR}
        $ENV{RDKIT_INCLUDE_PATH}
        $ENV{RDKIT_BASE}\\Code
        $ENV{RDBASE}\\Code
        C:\\RDKit\\include
        C:\\RDKit\\Code
      )
      find_path(RDKIT_INCLUDE_EXT_DIR INCHI-API/inchi.h
        PATHS
        ${RDKIT_DIR}\\External
        $ENV{RDKIT_INCLUDE_EXT_DIR}
        $ENV{RDKIT_INCLUDE_EXT_PATH}
        $ENV{RDKIT_BASE}\\External
        $ENV{RDBASE}\\External
        C:\\RDKit\\include
        C:\\RDKit\\External
      )
    else()
      find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
        PATHS
          ${RDKIT_DIR}/Code
          $ENV{RDKIT_INCLUDE_DIR}
          $ENV{RDKIT_INCLUDE_PATH}
          $ENV{RDKIT_BASE}/Code
          $ENV{RDBASE}/Code
          /usr/local/rdkit/include/Code
          /usr/local/rdkit/include
          /usr/local/rdkit/Code
          ~/rdkit/Code
      )
      find_path(RDKIT_INCLUDE_EXT_DIR INCHI-API/inchi.h
        PATHS
          ${RDKIT_DIR}/External
          $ENV{RDKIT_INCLUDE_EXT_DIR}
          $ENV{RDKIT_INCLUDE_EXT_PATH}
          $ENV{RDKIT_BASE}/External
          $ENV{RDBASE}/External
          /usr/local/rdkit/include/External
          /usr/local/rdkit/include
          /usr/local/rdkit/External
          ~/rdkit/External
      )      
    endif()
    if(RDKIT_INCLUDE_DIR)
       message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_DIR}")
    endif()
    if(RDKIT_INCLUDE_EXT_DIR)
       message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_EXT_DIR}")
    endif()    
  endif()

  if(NOT RDKIT_LIBRARIES)
    find_library(FILEPARSERS_LIB NAMES FileParsers
      PATHS
        ${RDKIT_DIR}/lib
        $ENV{RDKIT_LIB_DIR}
        $ENV{RDKIT_LIB_PATH}
        $ENV{RDKIT_LIBRARIES}
        $ENV{RDKIT_BASE}/lib
        $ENV{RDBASE}/lib
        /usr/local/rdkit/lib
        ~/rdkit/lib
        $ENV{LD_LIBRARY_PATH}
    )
    if(FILEPARSERS_LIB)
       GET_FILENAME_COMPONENT(RDKIT_LIBRARY_DIR ${FILEPARSERS_LIB} PATH)
       message(STATUS "Found RDKit libraries at ${RDKIT_LIBRARY_DIR}")

      # Note that the order of the following libraries is significant!!
      find_library(SMILESPARSE_LIB NAMES SmilesParse
                                   HINTS ${RDKIT_LIBRARY_DIR})
      find_library(DEPICTOR_LIB NAMES Depictor
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(CHEMTRANS_LIB NAMES ChemTransforms
                                HINTS ${RDKIT_LIBRARY_DIR})								
      find_library(GRAPHMOL_LIB NAMES GraphMol
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(RDGEOMETRYLIB_LIB NAMES RDGeometryLib
                                     HINTS ${RDKIT_LIBRARY_DIR})
      find_library(RDGENERAL_LIB NAMES RDGeneral
                                 HINTS ${RDKIT_LIBRARY_DIR})
      find_library(SUBSTRUCT_LIB NAMES SubstructMatch
                                 HINTS ${RDKIT_LIBRARY_DIR})      
      find_library(GASTEIGER_LIB NAMES PartialCharges
                                 HINTS ${RDKIT_LIBRARY_DIR})  
      find_library(DATASTRUCT_LIB NAMES DataStructs
                                 HINTS ${RDKIT_LIBRARY_DIR}) 
      find_library(SUBGRAPH_LIB NAMES Subgraphs                                 
                                 HINTS ${RDKIT_LIBRARY_DIR})
      find_library(FINGERPRINT_LIB NAMES Fingerprints
                                 HINTS ${RDKIT_LIBRARY_DIR})  
      find_library(INCHI_LIB NAMES Inchi
                                 HINTS ${RDKIT_LIBRARY_DIR})      
      find_library(RDINCHI_LIB NAMES RDInchiLib
                                 HINTS ${RDKIT_LIBRARY_DIR}) 
	  find_library(OPT NAMES Optimizer
                                 HINTS ${RDKIT_LIBRARY_DIR}) 								 
	  find_library(FF NAMES ForceField
                                 HINTS ${RDKIT_LIBRARY_DIR})  
	  find_library(FFHELP NAMES ForceFieldHelpers
                                 HINTS ${RDKIT_LIBRARY_DIR}) 
	  find_library(CATALOG NAMES Catalogs
                                 HINTS ${RDKIT_LIBRARY_DIR})  
	  find_library(FRAGCAT NAMES FragCatalog
                                 HINTS ${RDKIT_LIBRARY_DIR})   

                                 
      set (RDKIT_LIBRARIES ${FILEPARSERS_LIB} ${SMILESPARSE_LIB}
              ${DEPICTOR_LIB} ${CHEMTRANS_LIB} ${GRAPHMOL_LIB} ${RDGEOMETRYLIB_LIB}
              ${RDGENERAL_LIB} ${SUBSTRUCT_LIB} ${GASTEIGER_LIB} 
              ${DATASTRUCT_LIB} ${SUBGRAPH_LIB} ${FINGERPRINT_LIB} 
              ${INCHI_LIB} ${RDINCHI_LIB} ${OPT} ${FF} ${FFHELP} ${CATALOG} ${FRAGCAT})
    endif()
    if(RDKIT_LIBRARIES)
            message(STATUS "Found RDKit library files at ${RDKIT_LIBRARIES}")
    endif()
  endif()

  if(RDKIT_INCLUDE_DIR AND RDKIT_INCLUDE_EXT_DIR AND RDKIT_LIBRARIES)
    set(RDKIT_FOUND TRUE)
  endif()

  mark_as_advanced(RDINCHI_LIB INCHI_LIB GASTEIGER_LIB SUBSTRUCT_LIB RDGENERAL_LIB RDGEOMETRYLIB_LIB GRAPHMOL_LIB DEPICTOR_LIB SMILESPARSE_LIB FILEPARSERS_LIB)
  mark_as_advanced(RDKIT_INCLUDE_DIR RDKIT_INCLUDE_EXT_DIR RDKIT_LIBRARIES)
endif()
