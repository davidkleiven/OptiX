IDIR=inc
ODIR=obj
SDIR=src
CXX=g++

#MEEP_IDIR = /home/david/Documents/MEEP/meep-1.3/src
#MEEP_LDIR = /home/david/Documents/MEEP/meep-1.3/install/lib
include paths.sh
HDF5_LDIR = /usr/lib/x86_64-linux-gnu/hdf5/serial

PLANE_WAVE := planeReflection.cpp dielectricSlab.cpp anisMat.cpp
PLANE_WAVE_OBJ := $(PLANE_WAVE:.cpp=.o)
PLANE_WAVE := $(addprefix $(SDIR)/, $(PLANE_WAVE))
PLANE_WAVE_OBJ := $(addprefix $(ODIR)/, $(PLANE_WAVE_OBJ))

STATIC_LIBS := -lhdf5 -lz -lgsl -lharminv -llapack -lcblas -latlas -lfftw3 -lm -ljsoncpp
SHARED_LIBS := -lmeep

planeReflection.out: $(PLANE_WAVE_OBJ)
	$(CXX) -o $@ $^ -L $(MEEP_LDIR) -L $(HDF5_LDIR) $(STATIC_LIBS) $(SHARED_LIBS)

subtractBkg.out: ${ODIR}/subtractBkg.o ${ODIR}/readCSVdata.o
	$(CXX) -o $@ $^

normalizeDFTFlux.out: ${ODIR}/normalizeDFTFlux.o ${ODIR}/readCSVdata.o
	$(CXX) -o $@ $^ -ljsoncpp

fourierPulse.out: ${ODIR}/fourierPulse.o ${ODIR}/readCSVdata.o
	$(CXX) -o $@ $^ -llapack -lcblas -lgsl -lm -ljsoncpp

transmittanceAngle.out: ${ODIR}/transmittanceAngle.o ${ODIR}/readCSVdata.o
	$(CXX) -o $@ $^

visualizePlane.out: obj/visualizePlane.o obj/dielectricSlab.o
	$(CXX) -o $@ $^ -L $(MEEP_LDIR) -L $(HDF5_LDIR) $(STATIC_LIBS) $(SHARED_LIBS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) -MMD -c -fPIC $< -o $@ -I $(IDIR) -I $(MEEP_IDIR)

clean:
	rm ${ODIR}/*.o

cleanExec:
	rm *.out
