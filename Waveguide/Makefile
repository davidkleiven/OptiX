IDIR=inc/
IDIR_HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include
LDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/lib
ODIR=obj
SDIR=src

SOURCES:= main.cpp cladding.cpp waveGuide.cpp waveGuideRadiusCurvature.cpp numerovSolver.cpp solver1D.cpp waveGuideStraight.cpp \
	        stdFDsolver.cpp crancNicholson.cpp waveGuideFDSimulation.cpp solver2D.cpp
OBJECTS := ${SOURCES:.cpp=.o}
SOURCES := $(addprefix ${SDIR}/, ${SOURCES})
OBJECTS := $(addprefix ${ODIR}/, ${OBJECTS})

TEST_SOURCES := thomasAlgTest.cpp numerovTest.cpp waveGuideRadiusCurvatureTest.cpp alltests.cpp testMatrixBuild.cpp
TEST_OBJECTS := ${TEST_SOURCES:%.cpp=%.o}
TEST_SOURCES := $(addprefix unittest/, ${TEST_SOURCES})
TEST_OBJECTS := $(addprefix unittest/, ${TEST_OBJECTS})

waveguide.out: ${OBJECTS}
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -lcblas -llapack -lgsl -ljsoncpp -lz -lhdf5_serial -lhdf5_serial_hl -lhdf5_cpp -lhdf5_hl_cpp

unittest.out: ${TEST_OBJECTS} $(filter-out obj/main.o, ${OBJECTS})
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -lcblas -lgsl -llapack -ljsoncpp -lz -lhdf5_serial -lhdf5_serial_hl -lhdf5_cpp -lhdf5_hl_cpp \
	-lboost_program_options -lboost_unit_test_framework

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -I ${IDIR_HDF5} -o $@

unittest/%.o: unittest/%.cpp
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -I ${IDIR_HDF5} -o $@

obj/alltest.o:
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -O ${SDIR} -I ${IDIR_HDF5} -o $@

clean:
	rm ${ODIR}/*.o *.out unittest/*.o
