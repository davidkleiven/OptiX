cmake_minimum_required (VERSION 2.6)

# Add libraries
include_directories(inc)
include_directories(/usr/include/jsoncpp)
add_subdirectory(src)
include_directories(/usr/local/lib)

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

find_package( HDF5 REQUIRED )
include_directories(${HDF5_INCLUDE_DIR})
set(LIB ${LIB} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

find_package( Armadillo REQUIRED )
include_directories(${Armadillo_INCLUDE_DIR})
set(LIB ${LIB} ${Armadillo_LIBRARIES})

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

# Add executables
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

add_executable( farFieldConvergence.out Applications/farFieldConvergence.cpp )
target_link_libraries( farFieldConvergence.out ${LIB} src )

add_executable( fresnelProp.out Applications/fresnelProp.cpp )
target_link_libraries( fresnelProp.out ${LIB} src )

add_executable( straightGuideSweepAngle.out Applications/straightGuideSweepAngle.cpp )
target_link_libraries( straightGuideSweepAngle.out ${LIB} src )

add_executable( testSaveArmaMat.out unittest/testSaveArmamat.cpp )
target_link_libraries( testSaveArmaMat.out ${LIB} src )

# Compile unittests
find_package( Boost COMPONENTS system filesystem unit_test_framework REQUIRED )
#find_package( Boost REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} )
add_executable( alltests.out unittest/alltests.cpp unittest/curvedWaveGuide2DTest.cpp unittest/testMatrixBuild.cpp
unittest/thomasAlgTest.cpp unittest/waveGuideRadiusCurvatureTest.cpp )
target_link_libraries( alltests.out ${LIB} src ${Boost_FILESYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
