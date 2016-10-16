#!/usr/bin/bash
# Macro for running the python plotFDPattern.py script for several UIDs

UIDFile="$1"
FNAME="data/singleCurvedWG"

while read -r line
do
  echo "Current UID: ${line}"
  python plotFDPattern.py --file="${FNAME}${line}.json"
done < "${UIDFile}"
