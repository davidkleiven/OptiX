import straightEigenmodes as se
import sys
import numpy as np
import h5py as h5

def main( argv ):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print ("Usage: python straightEigenModePerturbation.py --file=<h5file>")
            return
        else:
            print ("Unknown argument %s"%(argv))
            return

    if ( fname == "" ):
        print ("No filename given!")
        return

    eigModes = se.StraightEigenmodes()
    wavelength = 0.1569 # Should be included in the h5file
    #wavelength = 0.124 # Should be included in the h5file
    k0 = 2.0*np.pi/wavelength
    try:
        with h5.File( fname, 'r' ) as hf:
            eigModes.read( hf )
    except Exception as exc:
        print ("Problem when reading file!")
        print (str(exc))
        return

    eigModes.nPropagatingModes = 70
    eigModes.plotCorrection3Lowest()

if __name__ == "__main__":
    main(sys.argv[1:])
