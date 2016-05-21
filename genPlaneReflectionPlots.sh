# This script processes the output from the plane reflection simulation
DDIR="dataPlane"
FDIR="Figures/planeReflection"

mkdir -p ${FDIR}
echo "Plotting raw fields at monitor..."
# Generate plot for the case with no reflections
OFILE="${FDIR}/noBoundary.png"
IFILE="${DDIR}/NoBoundary/ezMonitor.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt

# Normal incidence
OFILE="${FDIR}/normalInc.png"
IFILE="${DDIR}/NormalInc/ezMonitor.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt
