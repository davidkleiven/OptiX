import sys
sys.path.append("../FresnelFDTD")
import matplotlib as mpl
import mplLaTeX as ml
mpl.rcParams.update(ml.params)
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colors 

MSG = "Usage: python gaussianBeam.py --lambda=<wavelength> --waist=<beamWaist> [--help]\n"
MSG += "help: Print this message\n"
MSG += "lamgda: wavelength\n"
MSG += "waist: Waist at focus\n"

def waist( z, w0, zR ):
    return w0*np.sqrt( 1.0 + (z/zR)**2 )

def radiusOfCurvature( z, zR ):
    return z*(1.0 + (zR/z)**2)

def inverseRadiusOfCurvature(z,zR):
    return z/(z**2 + zR**2)
 
def guoyPhase(z,zR):
    return np.arctan(z/zR)

def rayleighRange( w0, wavelength ):
    return np.pi*w0**2 /wavelength

def gaussianBeam( r, z, wavelength, beamWaist ):
    zR = rayleighRange(beamWaist, wavelength)
    k = 2.0*np.pi/wavelength
    invR = inverseRadiusOfCurvature(z,zR)
    w0 = waist(z,beamWaist,zR)
    return (beamWaist/w0)*np.exp(-(r/w0)**2)*np.exp(1j*(k*z+0.5*k*invR*r**2-guoyPhase(z,zR)))

def main(argv):
    N_WAVELENGTHS = 40
    fname = "Figures/gaussianBeam.png"
    for arg in argv:
        if ( arg.find("--help") != -1 ):
            print MSG
            return 0
        elif ( arg.find("--lambda=") != -1 ):
            wavelength = float(arg.split("--lambda=")[1])
        elif ( arg.find("--waist=") != -1 ):
            beamWaist = float(arg.split("--waist=")[1])
        else:
            print ("Unknown argument %s"%(arg))
            return 0            

    NPOINTS = 1001
    z = np.linspace(-0.5*N_WAVELENGTHS*wavelength, 0.5*N_WAVELENGTHS*wavelength, NPOINTS)
    wMax = np.max(waist(z,beamWaist,rayleighRange(beamWaist,wavelength)))
    x = np.linspace(-2.0*wMax, 2.0*wMax, NPOINTS)
    Z, X = np.meshgrid(z,x)
    beam = gaussianBeam(X,Z,wavelength,beamWaist).real
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.contourf(Z,X,beam, 200, cmap="coolwarm")
    ax.set_xlabel("$z (\lambda)$")
    ax.set_ylabel("$r (\lambda)$")
    fig.savefig(fname, bbox_inches="tight")
    
    print ("Figure written to %s"%(fname))
    
if __name__ == "__main__":
    main(sys.argv[1:])
