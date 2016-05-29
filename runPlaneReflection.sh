#!/usr/bin/env bash
# This script runs the plane reflection program with different parameters
RUN_ALL=true # If false only non-exsisting simulations will be started

BKGDIR="dataPlane/NoBoundary"
EPS_HIGH=2.25
DELTA_ANGLE=5
ANGLE_MAX=90
ANGLE=5

POLARIZATION="s"
 
DDIR=()
while [ ${ANGLE} -lt ${ANGLE_MAX} ];
do
  DIR="dataPlane/Inc${ANGLE}"
  if [ ${RUN_ALL} = true ]; then
    # Clean directory
    rm -f ${DIR}/*.h5 ${DIR}/*.png ${DIR}/*.csv ${DIR}/*.gif
    rm -f ${DIR}/bkg/*.h5 ${DIR}/bkg/*.png ${DIR}/bkg/*.csv ${DIR}/bkg/*.gif
    rm -f ${DIR}/scatter/*.h5 ${DIR}/scatter/*.png ${DIR}/scatter/*.csv ${DIR}/scatter/*.gif
    mkdir -p ${DIR}
    mkdir -p ${DIR}/bkg
    mkdir -p ${DIR}/scatter
    
    # Run simulation
    ./planeReflection.out ${DIR}/bkg 1.0 ${ANGLE} ${POLARIZATION}
    ./planeReflection.out ${DIR}/scatter ${EPS_HIGH} ${ANGLE} ${POLARIZATION}
  elif [ ! -d "$DIR" ]; then
    # Directory does not exist
    mkdir ${DIR}
    mkdir ${DIR}/bkg
    mkdir ${DIR}/scatter
    # Run simulation
    ./planeReflection.out ${DIR}/bkg 1.0 ${ANGLE} ${POLARIZATION}
    ./planeReflection.out ${DIR}/scatter ${EPS_HIGH} ${ANGLE} ${POLARIZATION}
  fi

  DDIR+=(${DIR})
  ANGLE=$[${ANGLE}+${DELTA_ANGLE}]
done

# Normalize fluxes
for ((i=0;i<${#DDIR[@]};i++));
do
  ./normalizeDFTFlux.out "${DDIR[$i]}/scatter/transmittedFlux.csv" "${DDIR[$i]}/bkg/transmittedFlux.csv" 

done
