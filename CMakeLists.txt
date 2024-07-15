########################################################
# cmake file for building Marlin ACTSTracking Package
# @author Karol Krizka
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################

# project name
PROJECT( ACTSTracking )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version

FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDCore )

FIND_PACKAGE( ROOT REQUIRED)

FIND_PACKAGE( Acts REQUIRED COMPONENTS Core PluginJson PluginTGeo )

INCLUDE(GNUInstallDirs)

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### DATA ###########################################################
ADD_SUBDIRECTORY(data)
ADD_SUBDIRECTORY(example)

### LIBRARY #################################################################

# sources
CONFIGURE_FILE(src/config.h.in config.h)
ADD_LIBRARY(ACTSTracking SHARED
  ./src/ACTSProcBase.cxx
  ./src/ACTSTruthTrackingProc.cxx
  ./src/ACTSTruthCKFTrackingProc.cxx
  ./src/ACTSSeededCKFTrackingProc.cxx
  ./src/ACTSDuplicateRemoval.cxx
  ./src/TrackTruthProc.cxx
  ./src/GeometryIdSelector.cxx
  ./src/GeometryIdMappingTool.cxx
  ./src/Helpers.cxx
  )
TARGET_INCLUDE_DIRECTORIES(ACTSTracking
  PRIVATE
    ./ACTSTracking
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Marlin_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(ACTSTracking
  ${ROOT_LIBRARIES} DD4hep::DDCore
  ActsCore ActsPluginTGeo ActsPluginJson)

INSTALL( TARGETS ACTSTracking DESTINATION lib )
INSTALL( DIRECTORY ACTSTracking DESTINATION include FILES_MATCHING PATTERN "*.hxx" )


