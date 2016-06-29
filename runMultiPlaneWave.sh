TOP_FOLDER=dataPlane
SUBFOLDER_BASE="VarResMultInc"
DDIR_BASE="${TOP_FOLDER}/${SUBFOLDER_BASE}"
EPS_HIGH=2.25
ANGLES=(5 20 45 75 85)
REL_BAND_WIDTHS=(0.5 0.5 0.5 0.05 0.01)
NFREQ=(200 200 200 20 20)
RESOLUTION=(10 10 10 30 30)
CONTROL_FILENAME="${DDIR_BASE}.txt"
VERBOSE=false
POLARISATIONS=("s" "p")

# Control filename variables
IS_SIMULATED="simulated"
IS_NORMALIZED="normalized"
IS_FFT="fourier"

echo "Using control file ${CONTROL_FILENAME}..."

# Check if ddir base exists
if ( ls ${TOP_FOLDER} | grep ${SUBFOLDER_BASE} )
then
  echo "Folders starting with ${DDIR_BASE} already exists."
  echo "All contents in these will be deleted."
  echo "Do you want to continue? (yes/no)"
  ANSWER="no"
  read ANSWER
  if [ ![${ANSWER} == "yes"] ]
  then
    exit 1
  fi
fi

# NOTE: It is important that the background run and the actual run is executed succesively.
# The actual run relies on a temporary file created by the background run.
# This file is overwritten everytime planeReflection is called wth a directory name containing "bkg"
for pol in ${POLARISATIONS[@]}
do
  for ((i=0;i<${#ANGLES[@]};i++));
  do
    CTR_MSG="${pol}${ANGLES[$i]}Eps${EPS_HIGH}"
    DDIR="${DDIR_BASE}${ANGLES[$i]}${pol}"
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

      ./planeReflection.out ${BKGDIR} 1.0 ${ANGLES[$i]} ${pol} ${REL_BAND_WIDTHS[$i]} ${NFREQ[$i]} ${RESOLUTION[$i]}

      # Compute with reflection
      mkdir -p ${ODIR}

      ./planeReflection.out ${ODIR} ${EPS_HIGH} ${ANGLES[$i]} ${pol} ${REL_BAND_WIDTHS[$i]} ${NFREQ[$i]} ${RESOLUTION[$i]}
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
      ./normalizeDFTFlux.out ${ODIR}/transmittedFlux.json ${BKGDIR}/transmittedFlux.json
      echo "${CTR_MSG}${IS_NORMALIZED}" >> "${CONTROL_FILENAME}"
    fi

    if grep -Fxq "${CTR_MSG}${IS_FFT}" "${CONTROL_FILENAME}"
    then
      if ${VERBOSE}
      then
        echo "Has already FFT the fields..."
      fi
    else
      ./fourierPulse.out ${ODIR}/realField.json ${BKGDIR}/realField.json ${ANGLES[$i]} 
      echo "${CTR_MSG}${IS_FFT}" >> "${CONTROL_FILENAME}"
    fi
  done
done
