#!/usr/bin/env bash
# This script runs the plane reflection program with different parameters
RUN_ALL=true # If false only non-exsisting simulations will be started

DDIR=( "dataPlane/NoBoundary" "dataPlane/NormalInc")
EPS_HIGH=(1.0 2.25)
ANGLE=(0.0 0.0)

POLARIZATION="s"
COMP_FFT=false
for ((i=0;i<${#DDIR[@]};i++));
do
  DIR=${DDIR[$i]}
  if [ ${RUN_ALL} = true ]; then
    # Clean directory
    rm -f ${DIR}/*.h5 ${DIR}/*.png ${DIR}/*.csv ${DIR}/*.gif
    
    # Run simulation
    ./planeReflection.out ${DIR} ${EPS_HIGH[$i]} ${ANGLE[$i]} ${POLARIZATION}
    COMP_FFT=true
  elif [ ! -d "$DIR" ]; then
    # Directory does not exist
    mkdir ${DIR}
    # Run simulation
    ./planeReflection.out ${DIR} ${EPS_HIGH[$i]} ${ANGLE[$i]} ${POLARIZATION}
    COMP_FFT=true
  fi

  if [ ${COMP_FFT} == true ]; then
    ./fourierPulse.out ${DIR}/ezMonitorTrans.csv
  fi
  COMP_FFT=false
done

# Normalize fluxes
for ((i=1;i<${#DDIR[@]};i++));
do
  ./normalizeDFTFlux.out "${DDIR[$i]}/ezMonitorTransFourier.csv" "${DDIR[0]}/ezMonitorTransFourier.csv"
done
