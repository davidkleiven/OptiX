#!/usr/bin/env bash
# This script runs the plane reflection program with different parameters
RUN_ALL=true # If false only non-exsisting simulations will be started

BKGDIR="dataPlane/NoBoundary"
EPS_HIGH=2.25
DELTA_ANGLE=5
ANGLE_MAX=80
ANGLE=0

POLARIZATION="s"
COMP_FFT=false
if [ ${RUN_ALL} == true ] || [ ! -d ${BKGDIR}]; then
    rm -f ${BKGDIR}/*.h5 ${BKGDIR}/*.png ${BKGDIR}/*.csv ${BKGDIR}/*.gif
    mkdir -p ${BKGDIR}
  ./planeReflection.out ${BKGDIR} 1.0 0.0 ${POLARIZATION}
  ./fourierPulse.out ${BKGDIR}/ezMonitorTrans.csv
fi
  
DDIR=()
while [ ${ANGLE} -lt ${ANGLE_MAX} ];
do
  DIR="dataPlane/Inc${ANGLE}"
  if [ ${RUN_ALL} = true ]; then
    # Clean directory
    rm -f ${DIR}/*.h5 ${DIR}/*.png ${DIR}/*.csv ${DIR}/*.gif
    mkdir -p ${DIR}
    
    # Run simulation
    ./planeReflection.out ${DIR} ${EPS_HIGH} ${ANGLE} ${POLARIZATION}
    COMP_FFT=true
  elif [ ! -d "$DIR" ]; then
    # Directory does not exist
    mkdir ${DIR}
    # Run simulation
    ./planeReflection.out ${DIR} ${EPS_HIGH} ${ANGLE[$i]} ${POLARIZATION}
    COMP_FFT=true
  fi

  if [ ${COMP_FFT} == true ]; then
    ./fourierPulse.out ${DIR}/ezMonitorTrans.csv
  fi
  COMP_FFT=false
  DDIR+=(${DIR})
  ANGLE=$[${ANGLE}+${DELTA_ANGLE}]
done

# Normalize fluxes
for ((i=0;i<${#DDIR[@]};i++));
do
  ./normalizeDFTFlux.out "${DDIR[$i]}/ezMonitorTransFourier.csv" "${BKGDIR}/ezMonitorTransFourier.csv"
  #./normalizeDFTFlux.out "${DDIR[$i]}/transmittedFlux.csv" "${BKGDIR}/transmittedFlux.csv"
done
