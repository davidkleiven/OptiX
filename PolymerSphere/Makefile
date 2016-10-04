SDIR=./
ODIR=./
IDIR=./
COMMON_LIBS=-lscuff -lcblas -llapack -ljsoncpp

include libHeaders.sh
#SCUFF_HEADER=/usr/local/include/scuff-em
#BUFF_HEADER=/usr/local/include/buff-em

sphereScat.out: ${ODIR}/sphereScattering.o
	g++ $^ -o $@ -L ${LIB_PATHS} ${COMMON_LIBS}

sphereVIE.out: ${ODIR}/sphereScatVIE.o
	g++ $^ -o $@ -L ${LIB_PATHS} ${LIBS} -lbuff

${ODIR}/%.o: ${SDIR}/%.cpp
	g++ -c -std=c++11 $< -o $@ -I ${SCUFF_HEADER} -I ${IDIR} -I ${BUFF_HEADER}

clean:
	rm ${ODIR}/*.o
