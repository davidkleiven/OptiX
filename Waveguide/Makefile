IDIR=inc
ODIR=obj
SDIR=src

SOURCES :=cladding.cpp waveGuide.cpp waveGuideRadiusCurvature.cpp main.cpp numerovSolver.cpp solver1D.cpp
OBJECTS := ${SOURCES:.cpp=.o}
SOURCES := $(addprefix ${SDIR}/, ${SOURCES})
OBJECTS := $(addprefix ${ODIR}/, ${OBJECTS})

waveguide.out: ${OBJECTS}
	g++ -std=c++11 -o $@ $^ -lcblas -llapack -lgsl

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -std=c++11 -MMD -c $< -I ${IDIR} -o $@

clean:
	rm ${ODIR}/*.o *.out
