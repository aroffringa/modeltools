cmake_minimum_required(VERSION 2.8.8)

project(modeltools)

add_definitions(-DAOPROJECT)

if(PORTABLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb -Wvla -Wall -DNDEBUG -std=c++11")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb -Wvla -Wall -Wextra -DNDEBUG -march=native -std=c++11")
endif(PORTABLE)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

SET(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
include_directories(${CASACORE_INCLUDE_DIRS})

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIR})

find_library(FFTW3_LIB fftw3 REQUIRED)
find_library(FFTW3_THREADS_LIB fftw3_threads REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)
include_directories(${FFTW3_INCLUDE_DIR})

#Prevent accidentally finding old BoostConfig.cmake file from casapy
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS date_time system filesystem program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

find_library(PTHREAD_LIB pthread REQUIRED)
find_library(DL_LIB dl REQUIRED)
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
include_directories(${GSL_INCLUDE_DIR})
find_library(GSL_CBLAS_LIB NAMES gslcblas)
add_definitions(-DHAVE_GSL)

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

add_executable(bbs2model bbs2model.cpp nlplfitter.cpp model/model.cpp)
target_link_libraries(bbs2model ${GSL_LIB} ${GSL_CBLAS_LIB})

add_executable(editmodel editmodel.cpp fitsiochecker.cpp fitsreader.cpp nlplfitter.cpp polynomialfitter.cpp rmsynthesis.cpp spectrumft.cpp model/model.cpp deconvolution/spectralfitter.cpp)
target_link_libraries(editmodel ${Boost_SYSTEM_LIBRARY} ${CFITSIO_LIBRARY} ${CASACORE_LIBRARIES} ${GSL_LIB} ${GSL_CBLAS_LIB} ${FFTW3_LIB})

add_executable(cluster cluster.cpp nlplfitter.cpp model/model.cpp)
target_link_libraries(cluster ${GSL_LIB} ${GSL_CBLAS_LIB})

add_executable(render render.cpp image.cpp nlplfitter.cpp fftconvolver.cpp fftwmanager.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp model/model.cpp modelrenderer.cpp system.cpp wsclean/logger.cpp)
target_link_libraries(render ${Boost_DATE_TIME_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB} ${CFITSIO_LIBRARY} ${CASACORE_LIBRARIES} ${FFTW3_LIB} ${FFTW3_THREADS_LIB})

message(STATUS "Flags passed to C++ compiler: " ${CMAKE_CXX_FLAGS})

