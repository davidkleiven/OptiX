
cmake_minimum_required (VERSION 2.6)

project( Waveguide )
set( CMAKE_CXX_STANDARD 11 )

include_directories( include )
set( CMAKE_LIBRARY_PATH /usr/lib/x86_64-linux-gnu )

set( USE_MPI TRUE )
if ( USE_MPI )
  MESSAGE( STATUS "Compiling applications with mpic++" )
  set(CMAKE_CXX_COMPILER mpic++)
else()
  MESSAGE( STATUS "Compiling applications with g++")
  set(CMAKE_CXX_COMPILER g++)
endif()

add_definitions( -O2 )
if ( NOT USE_MPI )
  MESSAGE( STATUS "Compiling with openmp flag")
  set( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fopenmp")
endif()

###################### FIND STATIC LIBRARIES ###################################
# Static HDF5 library
if( USE_MPI )
  #set( HDF5_USE_STATIC_LIBRARIES TRUE )
  MESSAGE( STATUS "Using the openMPI version of HDF5" )
  find_package( HDF5 REQUIRED COMPONENTS C CXX HL )

  #set( HDF5_ROOT "/usr/lib64/openmpi/lib/")
  #find_package( HDF5 REQUIRED COMPONENTS C CXX HL )
  #MESSAGE( STATUS ${HDF5_LIBRARIES})
  #if ( NOT HDF5_IS_PARALLEL )
  #  MESSAGE( FATAL_ERROR "Did not find HDF5 version with parallel IO support")
#endif()
  include_directories( ${HDF5_INCLUDE_DIRS} )

  set(LIB ${LIB} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES})
endif()

# Meep library
if ( USE_MPI )
  MESSAGE( STATUS "Compiling with MPI version of MEEP" )
  find_library( LIBMEEP meep_mpi )
else()
  MESSAGE( STATUS "Compiling with single core version of MEEP" )
  find_library( LIBMEEP meep  PATH "/usr/lib/x86_64-linux-gnu" )
endif()
if ( NOT LIBMEEP )
  MESSAGE( FATAL_ERROR "Did not find the MEEP library!")
endif()

find_path( MEEPH "meep.hpp" )
if ( NOT MEEPH )
  MESSAGE( FATAL_ERROR "Did not find the MEEP header files!" )
endif()

# Harminv library
find_library( LIBHARM harminv )
if ( NOT LIBHARM )
  MESSAGE( FATAL_ERROR "Did not find the harminv library!" )
endif()

find_path( HARMH "harminv.h" )
if ( NOT HARMH )
  MESSAGE( FATAL_ERROR "Did not find the harminv header files!" )
endif()

# Add subdirectory and create a static library from the source files
add_subdirectory(src)

################# FIND SHARED LIBRARIES ########################################
find_package( OpenMP )
if ( OpenMP_found )
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

# Lapack
find_package( LAPACK REQUIRED )
set( LIB ${LIB} ${LAPACK_LIBRARIES})

# FFTW library
find_library( FFTWLIB fftw3 )
if ( NOT FFTWLIB )
  MESSAGE( FATAL_ERROR "Did not find the FFTW3 library!")
endif()
set( LIB ${LIB} ${FFTWLIB})
# Armadillo library
find_package( Armadillo REQUIRED )
include_directories(${ARMADILLO_INCLUDE_DIRS})
set(LIB ${LIB} ${ARMADILLO_LIBRARIES})

# HDF5 library
if( NOT USE_MPI )
  MESSAGE( STATUS "Using the serial version of HDF5")
  find_package( HDF5 REQUIRED COMPONENTS C CXX HL)
  include_directories( ${HDF5_INCLUDE_DIRS} )

  set(LIB ${LIB} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES})
  MESSAGE( STATUS ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} )
endif()

# Find pthreads library
find_package( OpenMP )
if ( OpenMP_found )
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

find_package( OpenGL REQUIRED )


# Check for the VISA library
find_library( VISA_LIB visa HINTS /usr/local/lib )
if ( VISA_LIB )
  set( LIB ${LIB} ${VISA_LIB} )
  set( HAVE_LIB_VISA 1 )

  find_package(SFML REQUIRED system window graphics network audio)
  include_directories(SFML_INCLUDE_DIR)
  set( LIB ${LIB} ${SFML_LIBRARIES} )
else()
  MESSAGE ( STATUS "Did not find the VISA library")
  set( HAVE_LIB_VISA 0 )
endif ()

find_path( VISA_HEADER visa/ )
if ( NOT VISA_HEADER )
  MESSAGE(FATAL_ERROR "Did not find the VISA header files!")
  set( HAVE_LIB_VISA 0 )
endif( NOT VISA_HEADER )

include_directories( FLTK_INCLUDE_DIR )
set( LIB ${LIB} ${FLTK_BASE_LIBRARY} )
set( LIB ${LIB} ${FLTK_GL_LIBRARY} )
set( LIB ${LIB} ${FLTK_FORMS_LIBRARY} )

find_package( GSL REQUIRED )
set( LIB ${LIB} ${GSL_LIBRARIES} )

find_library( JSONCPP_LIB jsoncpp )
if ( JSONCPP_LIB )
  set( LIB ${LIB} ${JSONCPP_LIB})
else()
  MESSAGE( FATAL_ERROR "Did not find the jsoncpp library!")
endif()
find_path( JSON_HEADER json)
if ( NOT JSON_HEADER )
  MESSAGE( FATAL_ERROR "Did not find the json headers!")
endif( NOT JSON_HEADER )


configure_file( "include/config.h.in" "include/config.h" )
set( LIB ${LIB} sourcefiles)
add_executable( coccolith.out coccolith.cpp )
target_link_libraries( coccolith.out ${LIB} )

add_executable(steadyCoccolith.out steadyCoccolith.cpp )
target_link_libraries(steadyCoccolith.out ${LIB} )

add_executable(cwSolveMWE.out cwSolveMWE.cpp )
target_link_libraries(cwSolveMWE.out ${LIB} )
