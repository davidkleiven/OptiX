SDIR=src
ODIR=obj
LIBS=-lscuff -lcblas -llapack -ljsoncpp
SCUFF_HEADER=/usr/local/include/scuff-em

fresnelBEM.out: ${ODIR}/fresnelBEM.o 
	g++ $^ -o $@ ${LIBS}

${ODIR}/fresnelBEM.o: ${SDIR}/fresnelBEM.cpp
	g++ -c -std=c++11 $< -o $@ -I ${SCUFF_HEADER}

clean:
	rm ${ODIR}/*.o
