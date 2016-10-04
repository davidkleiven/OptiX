SDIR=./
ODIR=./
IDIR=./
LIBS=-lscuff -lbuff -lcblas -llapack -ljsoncpp
include libHeaders.sh
#SCUFF_HEADER=/usr/local/include/scuff-em
#BUFF_HEADER=/usr/local/include/buff-em

sphereScat.out: ${ODIR}/sphereScattering.o
	g++ $^ -o $@ -L ${LIB_PATHS} ${LIBS}

sphereVIE.out: ${ODIR}/sphereScatVIE.o
	g++ $^ -o $@ -L ${LIB_PATHS} ${LIBS}

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -c -std=c++11 $< -o $@ -I ${SCUFF_HEADER} -I ${IDIR} -I ${BUFF_HEADER}

clean:
	rm ${ODIR}/*.o
