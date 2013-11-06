##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-code/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
set ( BASE_HEADERS FragmentGraph.h Features.h Config.h
FragmentGraphGenerator.h FragmentTreeNode.h Inference.h IPFP.h MILP.h
MolData.h Param.h Util.h Message.h Identifier.h Comparators.h)

set ( TRAIN_HEADERS Comms.h EM.h )
 
set ( BASE_SRC_FILES  FragmentGraph.cpp 
                      Features.cpp
                      Config.cpp
                      FragmentGraphGenerator.cpp
                      FragmentTreeNode.cpp
					  Identifier.cpp
                      Inference.cpp
                      IPFP.cpp
					  Message.cpp 
                      MILP.cpp
                      MolData.cpp
                      Param.cpp
                      Util.cpp
					  Comparators.cpp )
                     
set ( TRAIN_SRC_FILES Comms.cpp 
                      EM.cpp )
                      
if (INCLUDE_TRAIN)
    set ( SRC_FILES ${BASE_SRC_FILES} ${TRAIN_SRC_FILES} )
    set ( HEADER_FILES ${BASE_HEADERS} ${TRAIN_HEADERS} )
    set ( REQUIRED_LIBS ${RDKIT_LIBRARIES} ${Boost_LIBRARIES}
${LPSOLVE_LIBRARY} ${MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})    
else ()
    set ( SRC_FILES ${BASE_SRC_FILES} )
    set ( HEADER_FILES ${BASE_HEADERS} )
    set ( REQUIRED_LIBS ${RDKIT_LIBRARIES} ${Boost_LIBRARIES}
${LPSOLVE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})       
endif () 

add_library ( cfm-code STATIC ${SRC_FILES} ${HEADER_FILES} )
target_link_libraries ( cfm-code ${REQUIRED_LIBS} )

