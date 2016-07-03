LIBS=-lscuff -lcblas -llapack -ljsoncpp
SCUFF_HEADER=/usr/local/include/scuff-em

fresnelBEM.out: fresnelBEM.o
	g++ $^ -o $@ ${LIBS}

%.o: %.cpp
	g++ -c -std=c++11 $< -o $@ -I ${SCUFF_HEADER}
