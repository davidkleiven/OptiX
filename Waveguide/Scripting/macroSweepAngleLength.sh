#!/usr/bin/bash
LENGTHS=(3.0 2.9 2.8 2.7 2.6 2.5 2.4 2.3 2.2)
DFOLDER="lengthSweep"
mkdir -p ${DFOLDER}
i=0
for wglength in ${LENGTHS[@]}
do
  .././incidentAngleSweep.out --wglength=${wglength} --dfolder=${DFOLDER} --noUID
  python ../incidentAngleSweep.py --file=${DFOLDER}/angleSweep.h5 --noshow --figname="${DFOLDER}/frame${i}.png"
  $i=$i+1
done

#cd ${DFOLDER}
#avconv -r 10 -i frame%d.png anim.mp4
