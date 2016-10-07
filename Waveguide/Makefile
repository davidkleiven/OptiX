IDIR=inc/
IDIR_HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include
LDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/lib
ODIR=obj
SDIR=src

SOURCES :=cladding.cpp waveGuide.cpp waveGuideRadiusCurvature.cpp main.cpp numerovSolver.cpp solver1D.cpp
OBJECTS := ${SOURCES:.cpp=.o}
SOURCES := $(addprefix ${SDIR}/, ${SOURCES})
OBJECTS := $(addprefix ${ODIR}/, ${OBJECTS})

waveguide.out: ${OBJECTS}
	g++ -std=c++11 -o $@ $^ -L ${LDIR} -lcblas -llapack -lgsl -ljsoncpp -lz -lhdf5_serial -lhdf5_serial_hl -lhdf5_cpp -lhdf5_hl_cpp

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -I ${IDIR_HDF5} -o $@

clean:
	rm ${ODIR}/*.o *.out
