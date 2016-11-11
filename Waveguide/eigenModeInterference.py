import sys
import h5py as h5
import eigenModes as em
import numpy as np

def amp(x):
    return 1.0

def main( argv ):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print ("Usage: python eigenModeInterference.py --file=<h5file>")
            return
        else:
            print ("Unknown argument %s"%(argv))
            return

    if ( fname == "" ):
        print ("No filename given!")
        return

    eigModes = em.Eigenmodes()
    wavelength = 0.1569 # Should be included in the h5file
    k0 = 2.0*np.pi/wavelength
    try:
        with h5.File( fname, 'r' ) as hf:
            eigModes.read( hf )
    except Exception as exc:
        print ("Problem when reading file!")
        print (str(exc))
        return

    eigModes.nPropagatingModes = 15
    eigModes.exportTransFname = "data/transmissionExport.h5"
    absorb = eigModes.effectiveAbsorption()
    prop = eigModes.propagationConstants( k0 )
    coeff = eigModes.computeInitialCoefficient( amp )
    eigModes.contour( coeff, prop, absorb, k0, 400E3 )
    eigModes.plotAbsorption( coeff, absorb, k0, 400E3 )
    eigModes.transmissionByIntegrateOverWG( coeff, prop, absorb, k0, 400E3 )

if __name__ == "__main__":
    main( sys.argv[1:] )
