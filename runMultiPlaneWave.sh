DDIR_BASE=dataPlane/MultInc
EPS_HIGH=2.25
ANGLES=(20 45 85)


for ((i=0;i<${#ANGLES[@]};i++));
do
  # Compute background
  DDIR="${DDIR_BASE}${ANGLES[$i]}"
  ODIR=${DDIR}/WithEps
  BKGDIR=${DDIR}/bkg
  rm -f ${ODIR}/*
  rm -f ${BKGDIR}/*
  mkdir -p ${BKGDIR}

  ./planeReflection.out ${BKGDIR} 1.0 ${ANGLES[$i]} s

  # Compute with reflection
  mkdir -p ${ODIR}

  ./planeReflection.out ${ODIR} ${EPS_HIGH} ${ANGLES[$i]} s

  # Normalize
  ./normalizeDFTFlux.out ${ODIR}/transmittedFlux.csv ${BKGDIR}/transmittedFlux.csv
done
