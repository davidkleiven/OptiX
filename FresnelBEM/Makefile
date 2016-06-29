SDIR=src
ODIR=obj
LIBS=-lscuff -lcblas -llapack
SCUFF_HEADER=/usr/local/include/scuff-em

fresnelBEM.out: ${ODIR}/fresnelBEM.o 
	g++ $^ -o $@ ${LIBS}

${ODIR}/fresnelBEM.o: ${SDIR}/fresnelBEM.cpp
	g++ -c $< -o $@ -I ${SCUFF_HEADER}
