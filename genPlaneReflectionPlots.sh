# This script processes the output from the plane reflection simulation
DDIR="dataPlane"
FDIR="Figures/planeReflection"

mkdir -p ${FDIR}
echo "Plotting raw fields at monitor..."
# Generate plot for the case with no reflections
OFILE="${FDIR}/noBoundarySource.png"
IFILE="${DDIR}/NoBoundary/ezMonitorSource.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt
OFILE="${FDIR}/noBoundaryTrans.png"
IFILE="${DDIR}/NoBoundary/ezMonitorTrans.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt

# Normal incidence
OFILE="${FDIR}/normalIncSource.png"
IFILE="${DDIR}/NormalInc/ezMonitorSource.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt

OFILE="${FDIR}/normalIncTrans.png"
IFILE="${DDIR}/NormalInc/ezMonitorTrans.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt

# Sinc functions
OFILE="${FDIR}/sincSourceBkg.png"
IFILE="${DDIR}/MultInc30/bkg/ezMonitorSource.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt

OFILE="${FDIR}/sincSourceWithEps.png"
IFILE="${DDIR}/MultInc30/WithEps/ezMonitorSource.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt

OFILE="${FDIR}/sincSourceWithEpsTrans.png"
IFILE="${DDIR}/MultInc30/WithEps/ezMonitorTrans.csv"
gnuplot -e "ofile='${OFILE}';ifile='${IFILE}'" plotFieldComponent.plt
