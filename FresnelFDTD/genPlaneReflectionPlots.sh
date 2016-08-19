# This script processes the output from the plane reflection simulation
DDIR="dataPlane/Anisotropic" # Assumes that continues with <angle><polarisation>. Example: MultInc20s
source anglesAnis.sh
FDIR="Figures/Anisotropic"

mkdir -p ${FDIR}
# Plot flux spectrum
echo "Plotting flux spectrum..."
python plotFlux.py "${DDIR}" "${FDIR}" "${ANGLES[@]}"

# Plot amplitudes
echo "Plotting amplitudes..."
python plotFieldCoeff.py "${DDIR}" "${FDIR}" "${ANGLES[@]}"

# Plot reflection angle
echo "Plotting reflection angle..."
python computeReflectionAngle.py "${DDIR}" "${FDIR}" "${ANGLES[@]}"

# Plot pulse
echo "Plotting the pulse ..."
python pulse.py "${DDIR}20s" "${FDIR}"
