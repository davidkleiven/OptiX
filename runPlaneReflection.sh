#!/usr/bin/env bash
# This script runs the plane reflection program with different parameters
RUN_ALL=false # If false only non-exsisting simulations will be started

DDIR=( "dataPlane/NoBoundary" "dataPlane/NormalInc" "dataPlane/Inc5" "dataPlane/Inc10")
EPS_HIGH=(1.0 2.25 2.25 2.25)
ANGLE=(0.0 0.0 5.0 10.0)

POLARIZATION="s"
for ((i=0;i<${#DDIR[@]};i++));
do
  if [ ${RUN_ALL} ]; then
    DIR=${DDIR[$i]}
    # Clean directory
    rm -f ${DIR}/*.h5 ${DIR}/*.png ${DIR}/*.csv ${DIR}/*.gif
    
    # Run simulation
    ./planeReflection.out ${DIR} ${EPS_HIGH[$i]} ${ANGLE[$i]} ${POLARIZATION}

  elif [! -d ${DIR}]
    # Directory does not exist
    mkdir ${DIR}
    # Run simulation
    ./planeReflection.out ${DIR} ${EPS_HIGH[$i]} ${ANGLE[$i]} ${POLARIZATION}
  fi
done

# Normalize fluxes
for ((i=1;i<${#DDIR[@]};i++));
do
  ./normalizeDFTFlux.out "${DDIR[$i]}/transmittedFlux.csv" "${DDIR[0]}/transmittedFlux.csv"
done
