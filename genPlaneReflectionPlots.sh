# This script processes the output from the plane reflection simulation
DDIR="dataPlane"
FDIR="Figures"

# Plot flux spectrum
echo "Plotting flux spectrum..."
python plotFlux.py

# Plot amplitudes
echo "Plotting amplitudes..."
python plotFieldCoeff.py

# Plot reflection angle
echo "Plotting reflection angle..."
python computeReflectionAngle.py
