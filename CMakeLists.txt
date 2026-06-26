##############################################################################
# CMake file for Experiment
##############################################################################

cmake_minimum_required(VERSION 3.5)

project(Experiment)

set(EXPERIMENT_WITH_CGAL TRUE)

if(EXPERIMENT_WITH_CGAL)
   find_package(CGAL REQUIRED)
   find_package(Eigen3 REQUIRED)
   include_directories("/usr/include/eigen3")
   add_definitions(-DEXPERIMENT_WITH_CGAL)
endif()

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

if(EXPERIMENT_WITH_CGAL)
   aux_source_directories(SOURCES "Source Files\\CGAL" CGAL)
endif()

gomgen(Experiment)


add_library(Experiment ${SOURCES})

target_link_libraries(Experiment mesh_gfx CGAL)

if(EXPERIMENT_WITH_CGAL)
   target_link_libraries(Experiment CGAL)
endif()
