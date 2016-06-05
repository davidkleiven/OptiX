DDIR_BASE=dataPlane/MultInc
EPS_HIGH=2.25
ANGLES=(20 45 75 85)
REL_BAND_WIDTHS=(0.5 0.5 0.05 0.01)
NFREQ=(200 200 20 20)
CONTROL_FILENAME="dataPlane/control.txt"
VERBOSE=false

# Control filename variables
IS_SIMULATED="simulated"
IS_NORMALIZED="normalized"
IS_FFT="fourier"

# NOTE: It is important that the background run and the actual run is executed succesively.
# The actual run relies on a temporary file created by the background run.
# This file is overwritten everytime planeReflection is called wth a directory name containing "bkg"
for ((i=0;i<${#ANGLES[@]};i++));
do
  CTR_MSG="s${ANGLES[$i]}Eps${EPS_HIGH}"
  DDIR="${DDIR_BASE}${ANGLES[$i]}"
  ODIR=${DDIR}/WithEps
  BKGDIR=${DDIR}/bkg
  if grep -Fxq "${CTR_MSG}${IS_SIMULATED}" "${CONTROL_FILENAME}" 
  then
    if ${VERBOSE}
    then
      echo "Nothing to simulate for ${ANGLES[$i]} degrees s-polarization"
    fi
  else
    # Compute background
    rm -f ${ODIR}/*
    rm -f ${BKGDIR}/*
    mkdir -p ${BKGDIR}

    ./planeReflection.out ${BKGDIR} 1.0 ${ANGLES[$i]} s ${REL_BAND_WIDTHS[$i]} ${NFREQ[$i]}

    # Compute with reflection
    mkdir -p ${ODIR}

    ./planeReflection.out ${ODIR} ${EPS_HIGH} ${ANGLES[$i]} s ${REL_BAND_WIDTHS[$i]} ${NFREQ[$i]}
    echo "${CTR_MSG}${IS_SIMULATED}" >> "${CONTROL_FILENAME}"
  fi

  if grep -Fxq "${CTR_MSG}${IS_NORMALIZED}" "${CONTROL_FILENAME}" 
  then
    if ${VERBOSE}
    then
      echo "Already normalized..."
    fi
  else
    # Normalize
    ./normalizeDFTFlux.out ${ODIR}/transmittedFlux.csv ${BKGDIR}/transmittedFlux.csv
    echo "${CTR_MSG}${IS_NORMALIZED}" >> "${CONTROL_FILENAME}"
  fi

  if grep -Fxq "${CTR_MSG}${IS_FFT}" "${CONTROL_FILENAME}"
  then
    if ${VERBOSE}
    then
      echo "Has already FFT the fields..."
    fi
  else
    ./fourierPulse.out ${ODIR}/realField.csv ${BKGDIR}/realField.csv ${ANGLES[$i]} 
    echo "${CTR_MSG}${IS_FFT}" >> "${CONTROL_FILENAME}"
  fi
done

# Run p-polrazation
for ((i=0;i<${#ANGLES[@]};i++));
do
  CTR_MSG="p${ANGLES[$i]}Eps${EPS_HIGH}"
  if grep -Fxq "${CTR_MSG}${IS_SIMULATED}" "${CONTROL_FILENAME}" 
  then
    if ${VERBOSE}
    then
      echo "Nothing to simulate for ${ANGLES[$i]} degrees p-polarization"
    fi
  else
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
    echo "${CTR_MSG}${IS_SIMULATED}" >> "${CONTROL_FILENAME}"
  fi

  if grep -Fxq "${CTR_MSG}${IS_NORMALIZED}" "${CONTROL_FILENAME}" 
  then
    if ${VERBOSE}
    then
      echo "Already normalized..."
    fi
  else
  # Normalize
    ./normalizeDFTFlux.out ${ODIR}/transmittedFlux.csv ${BKGDIR}/transmittedFlux.csv
    echo "${CTR_MSG}${IS_NORMALIZED}" >> "${CONTROL_FILENAME}"
  fi

  if grep -Fxq "${CTR_MSG}${IS_FFT}" "${CONTROL_FILENAME}"
  then
    if ${VERBOSE}
    then
      echo "Has already FFT the fields..."
    fi
  else
    ./fourierPulse.out ${ODIR}/realField.csv ${BKGDIR}/realField.csv ${ANGLES[$i]} 
    echo "${CTR_MSG}${IS_FFT}" >> "${CONTROL_FILENAME}"
  fi
done
