project(parma)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
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
# if not found parma will be compiled for serial execution
# Uncomment the find of MPI to compile and run parma serially
###############################################################################

#find_package(MPI)
if(MPI_FOUND)
  message(STATUS "Package MPI found.")
else(MPI_FOUND)
  message(STATUS "The package MPI was not found.")
  set(CMAKE_CXX_FLAGS "-DSERIAL")
endif(MPI_FOUND)

###############################################################################
# Link and include
###############################################################################

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
