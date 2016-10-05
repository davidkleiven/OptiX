include searchPaths.sh

.PHONY: clean

sphereSubstrate.out: sphereSubstrate.o
	g++ -std=c++11 -o $@ $^ -L ${BA_LIB} -l:_libBornAgainCore.so -l:_libBornAgainFit.so
%.o: %.cpp
	g++ -std=c++11 -MMD -c $< -I ${BA_HEADERS} -o $@

clean:
	rm *.o *.out *.d
