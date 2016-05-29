DDIR_BASE=dataPlane/MultInc
EPS_HIGH=2.25
ANGLES=(20 45 75 85)
REL_BAND_WIDTHS=(0.5 0.5 0.05 0.01)
NFREQ=(200 200 200 200)


for ((i=0;i<${#ANGLES[@]};i++));
do
  # Compute background
  DDIR="${DDIR_BASE}${ANGLES[$i]}"
  ODIR=${DDIR}/WithEps
  BKGDIR=${DDIR}/bkg
  rm -f ${ODIR}/*
  rm -f ${BKGDIR}/*
  mkdir -p ${BKGDIR}

  ./planeReflection.out ${BKGDIR} 1.0 ${ANGLES[$i]} s ${REL_BAND_WIDTHS[$i]} ${NFREQ[$i]}

  # Compute with reflection
  mkdir -p ${ODIR}

  ./planeReflection.out ${ODIR} ${EPS_HIGH} ${ANGLES[$i]} s ${REL_BAND_WIDTHS[$i]} ${NFREQ[$i]}

  # Normalize
  ./normalizeDFTFlux.out ${ODIR}/transmittedFlux.csv ${BKGDIR}/transmittedFlux.csv
done

# Run p-polrazation
for ((i=0;i<${#ANGLES[@]};i++));
do
  # Compute background
  DDIR="${DDIR_BASE}${ANGLES[$i]}p"
  ODIR=${DDIR}/WithEps
  BKGDIR=${DDIR}/bkg
  rm -f ${ODIR}/*
  rm -f ${BKGDIR}/*
  mkdir -p ${BKGDIR}

  ./planeReflection.out ${BKGDIR} 1.0 ${ANGLES[$i]} p ${REL_BAND_WIDTHS[$i]} ${NFREQ[$i]}

  # Compute with reflection
  mkdir -p ${ODIR}

  ./planeReflection.out ${ODIR} ${EPS_HIGH} ${ANGLES[$i]} p ${REL_BAND_WIDTHS[$i]} ${NFREQ[$i]}

  # Normalize
  ./normalizeDFTFlux.out ${ODIR}/transmittedFlux.csv ${BKGDIR}/transmittedFlux.csv
done
