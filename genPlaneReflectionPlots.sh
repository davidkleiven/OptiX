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

echo "Subtracting off background..."
# Subtract of the background
BKGFILE="${DDIR}/NoBoundary/ezMonitor.csv"
./subtractBkg.out "${DDIR}/NormalInc/ezMonitor.csv" ${BKGFILE} "${DDIR}/NormalInc/ezMonitorDiff.csv"

# Plot DIFF files
echo "Plotting reflected fields..."
# Normal incidence
OFILE="${FDIR}/normalIncDiff.png"
IFILE="${DDIR}/NormalInc/ezMonitorDiff.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt
