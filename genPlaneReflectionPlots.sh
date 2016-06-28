# This script processes the output from the plane reflection simulation
DDIR="dataPlane/MultInc" # Assumes that continues with <angle><polarisation>. Example: MultInc20s
INC_ANGLE=(5 20 45 75 85)
FDIR="Figures"

# Plot flux spectrum
echo "Plotting flux spectrum..."
python plotFlux.py ${DDIR} ${FDIR} ${INC_ANGLE[@]}

# Plot amplitudes
echo "Plotting amplitudes..."
python plotFieldCoeff.py

# Plot reflection angle
echo "Plotting reflection angle..."
python computeReflectionAngle.py

# Plot pulse
echo "Plotting the pulse ..."
python pulse.py
