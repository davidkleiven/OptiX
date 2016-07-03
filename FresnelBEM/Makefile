SDIR=src
ODIR=obj
IDIR=inc
LIBS=-lscuff -lcblas -llapack -ljsoncpp
SCUFF_HEADER=/usr/local/include/scuff-em

fresnelBEM.out: ${ODIR}/fresnelBEM.o ${ODIR}/FluxIntegrator.o
	g++ $^ -o $@ ${LIBS}

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -c -std=c++11 $< -o $@ -I ${SCUFF_HEADER} -I ${IDIR}

clean:
	rm ${ODIR}/*.o
