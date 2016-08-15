
IDIR=inc
ODIR=obj
SDIR=src
CXX=g++

# There is a config script that creates the file paths.sh in the FresnelFDTD directory
include ../FresnelFDTD/paths.sh
HDF5_LDIR = /usr/lib/x86_64-linux-gnu/hdf5/serial

STATIC_LIBS := -lhdf5 -lz -lgsl -lharminv -llapack -lcblas -latlas -lfftw3 -lm -ljsoncpp
SHARED_LIBS := -lmeep

slabguide.out: obj/slabguide.o
	$(CXX) -o $@ $^ -L $(MEEP_LDIR) -L $(HDF5_LDIR) $(STATIC_LIBS) $(SHARED_LIBS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) -MMD -c -fPIC $< -o $@ -I $(IDIR) -I $(MEEP_IDIR)

clean:
	rm ${ODIR}/*.o

cleanExec:
	rm *.out
