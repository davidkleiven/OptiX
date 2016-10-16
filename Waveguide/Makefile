IDIR=inc/
IDIR_HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include
LDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/lib
ARMA_PATH=/home/david/Documents/localLib/armadillo-7.400.4/
ARMA_INC=/home/david/Documents/localLib/include
ODIR=obj
SDIR=src
OPT=-O2#-pg

SOURCES:= main.cpp cladding.cpp waveGuide.cpp waveGuideRadiusCurvature.cpp numerovSolver.cpp solver1D.cpp waveGuideStraight.cpp \
	        stdFDsolver.cpp crancNicholson.cpp waveGuideFDSimulation.cpp solver2D.cpp curvedWaveGuide2D.cpp controlFile.cpp \
					straightWG2D.cpp
OBJECTS := ${SOURCES:.cpp=.o}
SOURCES := $(addprefix ${SDIR}/, ${SOURCES})
OBJECTS := $(addprefix ${ODIR}/, ${OBJECTS})

TEST_SOURCES := thomasAlgTest.cpp numerovTest.cpp waveGuideRadiusCurvatureTest.cpp alltests.cpp testMatrixBuild.cpp \
                curvedWaveGuide2DTest.cpp
TEST_OBJECTS := ${TEST_SOURCES:%.cpp=%.o}
TEST_SOURCES := $(addprefix unittest/, ${TEST_SOURCES})
TEST_OBJECTS := $(addprefix unittest/, ${TEST_OBJECTS})

waveguide.out: ${OBJECTS}
	g++ -std=c++11 ${OPT} -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} -lcblas -llapack -lgsl -ljsoncpp -lz -lhdf5_serial -lhdf5_serial_hl \
	-lhdf5_cpp -lhdf5_hl_cpp -larmadillo

unittest.out: ${TEST_OBJECTS} $(filter-out obj/main.o, ${OBJECTS})
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} -lcblas -lgsl -llapack -ljsoncpp -lz -lhdf5_serial -lhdf5_serial_hl \
	-lhdf5_cpp -lhdf5_hl_cpp -lboost_program_options -lboost_unit_test_framework

curvedWG2D.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/curvedWGCartesian.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} -lcblas -lgsl -llapack -ljsoncpp -lz -lhdf5_serial -lhdf5_serial_hl -lhdf5_cpp \
	-lhdf5_hl_cpp -larmadillo

convergence.out: $(filter-out obj/main.o, ${OBJECTS}) ${ODIR}/convergence.o
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -L ${ARMA_PATH} -lcblas -lgsl -llapack -ljsoncpp -lz -lhdf5_serial -lhdf5_serial_hl -lhdf5_cpp \
	-lhdf5_hl_cpp -larmadillo

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -std=c++11 -MMD ${OPT} -c $< -I ${IDIR} -I ${IDIR_HDF5} -I ${ARMA_INC} -o $@

unittest/%.o: unittest/%.cpp
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -I ${IDIR_HDF5} -I ${ARMA_INC} -o $@

obj/alltest.o:
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -O ${SDIR} -I ${IDIR_HDF5} -o $@

clean:
	rm ${ODIR}/*.o *.out unittest/*.o
