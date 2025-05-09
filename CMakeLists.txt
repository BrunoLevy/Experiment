##############################################################################
# CMake file for Experiment
##############################################################################

cmake_minimum_required(VERSION 3.5)

project(Experiment)

if(EXISTS ${CMAKE_SOURCE_DIR}/../../../../Vorpaline/trunk/)
   set(
      GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../../../../Vorpaline/trunk"
      CACHE PATH "full path to the Geogram (or Vorpaline) installation"
   )
else()
   set(
      GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../../../../geogram/"
      CACHE PATH "full path to the Geogram (or Vorpaline) installation"
   )
endif()

set(
   GRAPHITE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../../GraphiteThree/"
   CACHE PATH "full path to the Graphite installation"
)

include(${GRAPHITE_SOURCE_DIR}/cmake/graphite.cmake)

if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   message(STATUS "Using options file: ${CMAKE_SOURCE_DIR}/CMakeOptions.txt")
   include(${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
endif()


##############################################################################

include_directories(${CMAKE_SOURCE_DIR}/../../)

# DLLs/shared objects are generated in Graphite binary directory

set(
   CMAKE_RUNTIME_OUTPUT_DIRECTORY
   ${GRAPHITE_SOURCE_DIR}/${RELATIVE_BIN_OUTPUT_DIR}
)

set(
   CMAKE_LIBRARY_OUTPUT_DIRECTORY
   ${GRAPHITE_SOURCE_DIR}/${RELATIVE_LIB_OUTPUT_DIR}
)

##############################################################################

aux_source_directories(SOURCES "Source Files\\common" common)
aux_source_directories(SOURCES "Source Files\\algo" algo)
aux_source_directories(SOURCES "Source Files\\tools" tools)
aux_source_directories(SOURCES "Source Files\\commands" commands)
aux_source_directories(SOURCES "Source Files\\third_party\\BasicInterval" third_party/BasicInterval)
gomgen(Experiment)

add_library(Experiment ${SOURCES})

target_link_libraries(Experiment mesh_gfx)
