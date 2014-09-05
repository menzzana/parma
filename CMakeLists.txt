project(parma)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST} global.h)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

###############################################################################
# Set version
###############################################################################

set(PARMA_VERSION_MAJOR "0")
set(PARMA_VERSION_MINOR "1")
set(PARMA_VERSION_PATCH "0")
configure_file( version_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/version_config.h )

###############################################################################
# Add mpi
###############################################################################

find_package(MPI)
if(MPI_FOUND)
  message(STATUS "Package MPI found.")
else(MPI_FOUND)
  message(STATUS "The package MPI was not found.")
endif(MPI_FOUND)

###############################################################################
# Link and include
###############################################################################

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${MPI_C_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES})
