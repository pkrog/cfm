##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# CMakeLists.txt
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.
#
# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
cmake_minimum_required(VERSION 2.6)

project (cfm)

set (CFM_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE STRING "Where do you want the output binaries to install to?" )

#It is recommended, but not strictly necessary, to build and run the tests
set( INCLUDE_TESTS OFF CACHE BOOL "Do you want to build the tests?" )

#Building the training code requires extra dependencies, e.g. MPI
set( INCLUDE_TRAIN OFF CACHE BOOL "Do you want to build the training code?" )

# BOOST
set (Boost_ADDITIONAL_VERSIONS "1.51" "1.51.0" "1.47" "1.47.0")
set (Boost_USE_STATIC_LIBS ON)
include (FindBoost)
find_package ( Boost 1.47 COMPONENTS filesystem system serialization REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )

# RDKit
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} )
include (FindRDKit)
if (RDKIT_FOUND)
    include_directories( ${RDKIT_INCLUDE_DIR} ${RDKIT_INCLUDE_EXT_DIR} )
    link_directories ( ${RDKIT_LIBRARY_DIR} )
else ()
    message (FATAL_ERROR "Could not find RDKit" )
endif ()

# LPSolve
set ( LPSOLVE_INCLUDE_DIR "" CACHE STRING "Where are the LPSolve headers?" )
set ( LPSOLVE_LIBRARY_DIR "" CACHE STRING "Where are the LPSolve libraries?" )
include_directories ( ${LPSOLVE_INCLUDE_DIR} )
link_directories ( ${LPSOLVE_LIBRARY_DIR} )
find_library(LPSOLVE_LIBRARY NAMES lpsolve55 PATHS ${LPSOLVE_LIBRARY_DIR})
if (NOT LPSOLVE_LIBRARY)
  MESSAGE(FATAL_ERROR "Cannot find LPSOLVE libraries!")
endif (NOT LPSOLVE_LIBRARY)
mark_as_advanced( LPSOLVE_LIBRARY )

find_package( Threads )

#  MPI
if ( INCLUDE_TRAIN OR INCLUDE_TESTS )
    include (FindMPI)
    if ( MPI_FOUND )
        include_directories ( ${MPI_INCLUDE_PATH} )
    else ()
        message (FATAL_ERROR "Could not find MPI" )
    endif ()
endif ()

#Common Library
add_subdirectory( cfm-code )

#Executable for generating fragmentation graphs
include_directories( ${CMAKE_MODULE_PATH}/cfm-code )
add_subdirectory( fraggraph-gen )

#Executable for running tests
if (INCLUDE_TESTS)
    add_subdirectory( cfm-test )
endif ()

if (INCLUDE_TRAIN)
    add_subdirectory( cfm-train )
endif ()      

#Spectrum Predicter
add_subdirectory( cfm-predict )

#Identifier
add_subdirectory( cfm-id )

#Annotater
add_subdirectory( cfm-annotate )