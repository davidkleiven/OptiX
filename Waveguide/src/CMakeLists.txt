set( SOURCES cladding.cpp waveGuide.cpp waveGuideRadiusCurvature.cpp numerovSolver.cpp solver1D.cpp waveGuideStraight.cpp
	        stdFDsolver.cpp crancNicholson.cpp waveGuideFDSimulation.cpp solver2D.cpp curvedWaveGuide2D.cpp controlFile.cpp
					straightWG2D.cpp linearMap1D.cpp coupledCurvedWG.cpp paraxialSource.cpp planeWave.cpp
					gaussianBeam.cpp paraxialEquation.cpp cylindricalParaxialEquation.cpp curvedWGCylCrd.cpp
					beamSplitter.cpp gaussianWG.cpp linearRampWG.cpp fresnelPropagator.cpp exitFieldSource.cpp
					incidentAngleSweep.cpp borderTracker.cpp paraxialSimulation.cpp h5Attribute.cpp
				  refractiveIndex.cpp incidentAngleLengthSweep.cpp cylinder2D.cpp postProcessMod.cpp farFieldPost.cpp
				  transmittivity.cpp curvedWGConfMap.cpp multipleCurvedWG.cpp
				fixedValuesSource.cpp bendOptimizer.cpp fftSolver2D.cpp solver.cpp solver3D.cpp fftSolver3D.cpp
			postProcessing.cpp)

add_library( src STATIC ${SOURCES} )
