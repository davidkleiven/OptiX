
cmake_minimum_required (VERSION 2.6)

project( Waveguide )
set( CMAKE_CXX_STANDARD 11 )

include_directories( include )
add_subdirectory(src)
set( CMAKE_LIBRARY_PATH /usr/lib/x86_64-linux-gnu )

# Meep library
find_library( LIBMEEP meep  PATH "/usr/lib/x86_64-linux-gnu" )
if ( NOT LIBMEEP )
  MESSAGE( FATAL_ERROR "Did not find the MEEP library!")
endif()
set( LIB ${LIBMEEP} )

find_path( MEEPH "meep.hpp" )
if ( NOT MEEPH )
  MESSAGE( FATAL_ERROR "Did not find the MEEP header files!" )
endif()

# Harminv library
find_library( LIBHARM harminv )
if ( NOT LIBHARM )
  MESSAGE( FATAL_ERROR "Did not find the harminv library!" )
endif()
set( LIB ${LIB} ${LIBHARM} )
find_path( HARMH "harminv.h" )
if ( NOT HARMH )
  MESSAGE( FATAL_ERROR "Did not find the harminv header files!" )
endif()

# Armadillo library
find_package( Armadillo REQUIRED )
include_directories(${ARMADILLO_INCLUDE_DIRS})
set(LIB ${LIB} ${ARMADILLO_LIBRARIES})

# HDF5 library
find_package( HDF5 REQUIRED COMPONENTS C CXX HL)
include_directories( ${HDF5_INCLUDE_DIRS} )

MESSAGE( STATUS ${HDF5_INCLUDE_DIRS} )
set(LIB ${LIB} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES})
# Find the PEI library
find_library( PEI_LIB pei )
set( LIB ${LIB} ${PEI_LIB} )
if ( NOT PEI_LIB )
  MESSAGE( FATAL_ERROR "Did not find the pei library!")
endif ( NOT PEI_LIB )
find_path( PEI_HEADER pei/ )
if ( NOT PEI_HEADER )
  MESSAGE( FATAL_ERROR "Did not find the pei headers!")
endif ( NOT PEI_HEADER )

# Find pthreads library
find_package( OpenMP )
if ( OpenMP_found )
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

set ( FLTK_SKIP_FLUID True )
set( CMAKE_LIBRARY_PATH "/usr/lib/x86_64-linux-gnu")
#find_package( FLTK REQUIRED )

# Manual find the FLTK libraries
find_library( FLTK_BASE_LIBRARY fltk )
if( NOT FLTK_BASE_LIBRARY )
  MESSAGE( FATAL_ERROR "Did not find the fltk library")
endif()
find_library( FLTK_GL_LIBRARY fltk_gl)
if ( NOT FLTK_GL_LIBRARY )
  MESSAGE( FATAL_ERROR "Did not find the fltk_gl library!")
endif()
find_library( FLTK_FORMS fltk_forms )
if ( NOT FLTK_FORMS )
  MESSAGE( FATAL_ERROR "Did not find the fltk_forms library!")
endif()
find_library( FLTK_IMAGE_LIBRARY fltk_images)
if ( NOT FLTK_IMAGE_LIBRARY)
  MESSAGE( FATAL_ERROR "Did not find the fltk_images library!")
endif()
find_package( OpenGL REQUIRED )

find_package(SFML REQUIRED system window graphics network audio)
include_directories(SFML_INCLUDE_DIR)
set( LIB ${LIB} ${SFML_LIBRARIES} )

# Check for the VISA library
find_library( VISA_LIB visa HINTS /usr/local/lib )
if ( NOT VISA_LIB )
  MESSAGE (FATAL_ERROR "Did not find the VISA library")
endif ( NOT VISA_LIB )
set( LIB ${LIB} ${VISA_LIB} )

find_path( VISA_HEADER visa/ )
if ( NOT VISA_HEADER )
  MESSAGE(FATAL_ERROR "Did not find the VISA header files!")
endif( NOT VISA_HEADER )

include_directories( FLTK_INCLUDE_DIR )
set( LIB ${LIB} ${FLTK_BASE_LIBRARY} )
set( LIB ${LIB} ${FLTK_GL_LIBRARY} )
set( LIB ${LIB} ${FLTK_FORMS_LIBRARY} )

set( LIB ${LIB} sourcefiles)
add_executable( coccolith.out coccolith.cpp )
target_link_libraries( coccolith.out ${LIB} )
