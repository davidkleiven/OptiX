cmake_minimum_required (VERSION 2.6)

project( Waveguide )
set( CMAKE_CXX_STANDARD 11 )
#add_definitions( -DARMA_NO_DEBUG )
add_definitions( -O2 )
set( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fopenmp")

# Add libraries
include_directories(inc)
include_directories(core/include)
include_directories(/usr/include/jsoncpp)
include_directories( /usr/local/include )
include_directories( /usr/local/include/PaxPro )
include_directories( /usr/lib/x86_64-linux-gnu/hdf5/serial/include )

link_directories( /usr/local/lib )
link_directories( /usr/lib/x86_64-linux-gnu )
link_directories( ${CMAKE_CURRENT_SOURCE_DIR}/core/src )
set( CMAKE_LIBRARY_PATH /usr/lib/x86_64-linux-gnu )

# Find required packages
find_package( LAPACK REQUIRED )
include_directories(${LAPACK_INCLUDE_DIR})
set(LIB ${LAPACK_LIBRARIES})

find_package( BLAS REQUIRED )
include_directories(${BLAS_INCLUDE_DIR})
set(LIB ${LIB} ${BLAS_LIBRARIES})

find_package( GSL REQUIRED )
include_directories(${GSL_INCLUDE_DIR})
set(LIB ${LIB} ${GSL_LIBRARIES})

set(CMAKE_INCLUDE_PATH "/usr/include/jsoncpp")
find_package( jsoncpp REQUIRED )
include_directories(${JSONCPP_INCLUDE_DIR})
set(LIB ${LIB} ${JSONCPP_LIBRARIES})

# TODO: Seem to be something wrong witht the package detection
find_library( JSONCPP_LIB jsoncpp )
set( LIB ${LIB} jsoncpp.so)
if ( NOT JSONCPP_LIB )
  MESSAGE( FATAL_ERROR "Did not find the jsoncpp library!")
endif( NOT JSONCPP_LIB)
find_path( JSON_HEADER json)
if ( NOT JSON_HEADER )
  MESSAGE( FATAL_ERROR "Did not find the json headers!")
endif( NOT JSON_HEADER )

find_package( ZLIB REQUIRED )
include_directories(${ZLIB_INCLUDE_DIR})
set(LIB ${LIB} ${ZLIB_LIBRARIES})

find_package( HDF5 REQUIRED COMPONENTS C CXX HL)
include_directories( ${HDF5_INCLUDE_DIRS} )

MESSAGE( STATUS ${HDF5_INCLUDE_DIRS} )
set(LIB ${LIB} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES})

#set(CMAKE_INCLUDE_PATH "/usr/include")
find_package( Armadillo REQUIRED )
include_directories(${ARMADILLO_INCLUDE_DIRS})
set(LIB ${LIB} ${ARMADILLO_LIBRARIES})

find_package( SFML REQUIRED )
include_directories(${SFML_INCLUDE_DIR})
set(LIB ${LIB} libsfml-window.so libsfml-graphics.so libsfml-system.so)

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

# Check for the harminv library
find_library( HARMINV_LIB harminv HINTS /usr/local/lib )
set( LIB ${LIB}  ${HARMINV_LIB} )
if ( NOT HARMINV_LIB )
  MESSAGE(FATAL_ERROR "Did not find the Harminv library!")
endif( NOT HARMINV_LIB )
find_path( HARMINV_HEADER harminv.h )
if ( NOT HARMINV_HEADER )
  MESSAGE(FATAL_ERROR "Did not find the harminv header file!")
endif( NOT HARMINV_HEADER )

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

include_directories( FLTK_INCLUDE_DIR )
set( LIB ${LIB} ${FLTK_BASE_LIBRARY} )
set( LIB ${LIB} ${FLTK_GL_LIBRARY} )
set( LIB ${LIB} ${FLTK_FORMS_LIBRARY} )

# Find the FFTW3 library with threads
find_library( FFTW_OMP_LIB fftw3_omp )
if ( NOT FFTW_OMP_LIB )
  MESSAGE( FATAL_ERROR "Could not find the FFTW3 openMP library!")
endif()
set( LIB ${LIB} ${FFTW_OMP_LIB})

# Find the FFTW3 library
find_library( FFTW_LIB fftw3 )
if ( NOT FFTW_LIB )
  MESSAGE( FATAL_ERROR "Could not find the FFTW library")
endif()
find_path( FFTW_INC fftw3.h)
if ( NOT FFTW_INC )
  MESSAGE( FATAL_ERROR "Could not find the FFTW3 headers!")
endif()
set( LIB ${LIB} ${FFTW_LIB} )
set( LIB ${LIB} "m") # Add C-math library
include_directories( FFTW_INC )

find_library( PAXPRO NAMES libpaxpro.a paxpro )
if ( NOT PAXPRO )
  MESSAGE( FATAL_ERROR "Could not find the PaxPro library!" )
endif()
include_directories( "usr/local/include/PaXPro")

set( LIB ${LIB} ${PAXPRO} ) # It should find this library
add_subdirectory(src)

#======================== COMPILE EXECUTABLES ==================================
add_executable( curvedWG.out Applications/curvedWGCartesian.cpp )
target_link_libraries( curvedWG.out ${LIB} src )

add_executable( beamSplitSim.out Applications/beamSplitSim.cpp )
target_link_libraries( beamSplitSim.out ${LIB} src)

add_executable( borderDemo.out Applications/borderDemo.cpp )
target_link_libraries( borderDemo.out ${LIB} src )

add_executable( convergence.out Applications/convergence.cpp )
target_link_libraries( convergence.out ${LIB} src )

add_executable( coupledWaveGuideSimulation.out Applications/coupledWaveGuideSimulation.cpp )
target_link_libraries( coupledWaveGuideSimulation.out ${LIB} src )

add_executable( cylinderScattering.out Applications/cylinderScattering.cpp )
target_link_libraries( cylinderScattering.out ${LIB} src )

add_executable( fresnelProp.out Applications/fresnelProp.cpp )
target_link_libraries( fresnelProp.out ${LIB} src )

add_executable( straightGuideSweepAngle.out Applications/straightGuideSweepAngle.cpp )
target_link_libraries( straightGuideSweepAngle.out ${LIB} src )

add_executable( testSaveArmaMat.out unittest/testSaveArmamat.cpp )
target_link_libraries( testSaveArmaMat.out ${LIB} src )

add_executable( multiWG.out Applications/multiWG.cpp )
target_link_libraries( multiWG.out ${LIB} src )

add_executable( optimizeBend.out Applications/optimizeBend.cpp )
target_link_libraries( optimizeBend.out ${LIB} src )

add_executable( sphereScat.out Applications/sphereScattering.cpp )
target_link_libraries( sphereScat.out ${LIB} src )

add_executable( wg3DSim.out Applications/waveguide3Dsim.cpp )
target_link_libraries( wg3DSim.out ${LIB} src )

# Compile unittests
find_package( Boost COMPONENTS system filesystem unit_test_framework )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} )

if ( BOOST_FOUND )
  add_executable( alltests.out unittest/alltests.cpp unittest/curvedWaveGuide2DTest.cpp unittest/testMatrixBuild.cpp
  unittest/thomasAlgTest.cpp unittest/waveGuideRadiusCurvatureTest.cpp )
  target_link_libraries( alltests.out ${LIB} src ${Boost_FILESYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
endif()
