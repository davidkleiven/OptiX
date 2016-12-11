IDIR=inc/
IDIR_HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include
LDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/lib
ARMA_PATH=/home/david/Documents/localLib/armadillo-7.400.4/
ARMA_INC=/home/david/Documents/localLib/include
SFML_PATH=/usr/lib/x86_64-linux-gnu
APPDIR=Applications
ODIR=obj
SDIR=src
OPT=-O2#-pg

LIBS:=-lcblas -llapack -lgsl -ljsoncpp -lz -lhdf5_serial -lhdf5_serial_hl \
      -lhdf5_cpp -lhdf5_hl_cpp -larmadillo -lharminv
GRAPHLIBS := -lsfml-graphics -lsfml-window -lsfml-system -lvisa

SOURCES:= main.cpp cladding.cpp waveGuide.cpp waveGuideRadiusCurvature.cpp numerovSolver.cpp solver1D.cpp waveGuideStraight.cpp \
	        stdFDsolver.cpp crancNicholson.cpp waveGuideFDSimulation.cpp solver2D.cpp curvedWaveGuide2D.cpp controlFile.cpp \
					straightWG2D.cpp linearMap1D.cpp coupledCurvedWG.cpp harmonicInversion.cpp paraxialSource.cpp planeWave.cpp \
					gaussianBeam.cpp paraxialEquation.cpp cylindricalParaxialEquation.cpp curvedWGCylCrd.cpp \
					beamSplitter.cpp gaussianWG.cpp linearRampWG.cpp fresnelPropagator.cpp exitFieldSource.cpp \
					incidentAngleSweep.cpp visualizer.cpp colormaps.cpp visualizer1D.cpp borderTracker.cpp \
					refractiveIndex.cpp incidentAngleLengthSweep.cpp

OBJECTS := ${SOURCES:.cpp=.o}
SOURCES := $(addprefix ${SDIR}/, ${SOURCES})
OBJECTS := $(addprefix ${ODIR}/, ${OBJECTS})

TEST_SOURCES := thomasAlgTest.cpp numerovTest.cpp waveGuideRadiusCurvatureTest.cpp alltests.cpp testMatrixBuild.cpp \
                curvedWaveGuide2DTest.cpp
TEST_OBJECTS := ${TEST_SOURCES:%.cpp=%.o}
TEST_SOURCES := $(addprefix unittest/, ${TEST_SOURCES})
TEST_OBJECTS := $(addprefix unittest/, ${TEST_OBJECTS})

waveguide.out: ${OBJECTS}
	g++ -std=c++11 ${OPT} -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS} ${GRAPHLIBS}

unittest.out: ${TEST_OBJECTS} $(filter-out obj/main.o, ${OBJECTS})
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS} -lboost_program_options -lboost_unit_test_framework

curvedWG2D.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/curvedWGCartesian.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS} ${GRAPHLIBS}

convergence.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/convergence.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS} ${GRAPHLIBS}

wgHarmInv.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/wgHarmonicInversion.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS}

eigenProject.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/projectOntoEigenmodes.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS}

coupledWG.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/coupledWaveGuideSimulation.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS}

beamSplitSim.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/beamSplitSim.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS} ${GRAPHLIBS}

fresnelProp.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/fresnelProp.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS} ${GRAPHLIBS}

incidentAngleSweep.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/straightGuideSweepAngle.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} ${LIBS} ${GRAPHLIBS}

testVisualizer.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/testVisualizer.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} -L ${SFML_PATH} ${LIBS} ${GRAPHLIBS}

farFieldConvergence.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/farFieldConvergence.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} -L ${SFML_PATH} ${LIBS} ${GRAPHLIBS}

borderDemo.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/borderDemo.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} -L ${SFML_PATH} ${LIBS} ${GRAPHLIBS}

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -std=c++11 -MMD ${OPT} -c $< -I ${IDIR} -I ${IDIR_HDF5} -I ${ARMA_INC} -o $@

${ODIR}/%.o: ${APPDIR}/%.cpp
	g++ -std=c++11 -MMD ${OPT} -c $< -I ${IDIR} -I ${IDIR_HDF5} -I ${ARMA_INC} -o $@

unittest/%.o: unittest/%.cpp
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -I ${IDIR_HDF5} -I ${ARMA_INC} -o $@

obj/alltest.o:
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -O ${SDIR} -I ${IDIR_HDF5} -o $@

clean:
	rm ${ODIR}/*.o *.out unittest/*.o
