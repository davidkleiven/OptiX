#!/usr/bin/bash
LENGTHS=(3.0 2.95 2.9 2.85 2.8 2.75 2.7 2.65 2.6 2.55 2.5 2.45 2.4 2.35 2.3 2.25 2.2 2.15 2.1 2.05 2.0)
DFOLDER="lengthSweep"
mkdir -p ${DFOLDER}
i=0
for wglength in ${LENGTHS[@]}
do
  ./incidentAngleSweep.out --wglength=${wglength} --dfolder=${DFOLDER} --noUID --visualize
  python incAngleSweep.py --file="${DFOLDER}/angleSweep.h5" --noshow --figname="${DFOLDER}/frame${i}.png" --title="L=$wglength"
  i=$((i+1))
done

#cd ${DFOLDER}
#avconv -r 10 -i frame%d.png anim.mp4
