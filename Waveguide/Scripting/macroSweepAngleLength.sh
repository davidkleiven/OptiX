#!/usr/bin/bash
LENGTHS=(1.15 1.14 1.13 1.12 1.11 1.1 1.09 1.08 1.07 1.06 1.05 1.04 1.03 1.02 1.0 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.9)
DFOLDER="lengthSweep"
mkdir -p ${DFOLDER}
i=0
for wglength in ${LENGTHS[@]}
do
  ./incidentAngleSweep.out --wglength=${wglength} --dfolder=${DFOLDER} --noUID --visualize
  python incAngleSweep.py --file="${DFOLDER}/angleSweep.h5" --noshow --figname="${DFOLDER}/frame${i}.png" --title="L=$wglength"
  i=$((i+1))
done

DFOLDER="lengthSweepAlc"
mkdir -o ${DFOLDER}
for wglength in ${LENGTHS[@]}
do
  ./incidentAngleSweep.out --wglength=${wglength} --dfolder=${DFOLDER} --noUID --visualize --usealc
  python incAngleSweep.py --file="${DFOLDER}/angleSweep.h5" --noshow --figname="${DFOLDER}/frame${i}.png" --title="L=$wglength"
  i=$((i+1))
done

#cd ${DFOLDER}
#avconv -r 10 -i frame%d.png anim.mp4
