SDIR=./
ODIR=./
IDIR=./
LIBS=-lscuff -lcblas -llapack -ljsoncpp
SCUFF_HEADER=/usr/local/include/scuff-em

sphereScat.out: ${ODIR}/sphereScattering.o
	g++ $^ -o $@ ${LIBS}

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -c -std=c++11 $< -o $@ -I ${SCUFF_HEADER} -I ${IDIR}

clean:
	rm ${ODIR}/*.o
