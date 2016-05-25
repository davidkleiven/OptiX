DDIR=dataPlane/MultInc30
EPS_HIGH=2.25
ANGLE=30

# Compute background
BKGDIR=${DDIR}/bkg
mkdir -p ${BKGDIR}
rm -f ${BKGDIR}/*

./planeReflection.out ${BKGDIR} 1.0 ${ANGLE} s

# Compute with reflection
ODIR=${DDIR}/WithEps
mkdir -p ${ODIR}
rm -f ${ODIR}/*

./planeReflection.out ${ODIR} ${EPS_HIGH} ${ANGLE} s

# Normalize
./normalizeDFTFlux.out ${ODIR}/transmittedFlux.csv ${BKGDIR}/transmittedFlux.csv
