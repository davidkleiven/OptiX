set( SOURCES beamSplitter.cpp bendOptimizer.cpp coupledCurvedWG.cpp curvedWGConfMap.cpp
	curvedWGCylCrd.cpp cylinder2D.cpp cylindricalParaxialEquation.cpp fresnelPropagator.cpp
	gaussianWG.cpp incidentAngleLengthSweep.cpp
	incidentAngleSweep.cpp linearRampWG.cpp main.cpp multipleCurvedWG.cpp
	numerovSolver.cpp sphere.cpp straightWG2D.cpp
	waveguide3DPaths.cpp waveGuideRadiusCurvature.cpp waveguideShapes.cpp
waveGuideStraight.cpp )

add_library( src STATIC ${SOURCES} )
