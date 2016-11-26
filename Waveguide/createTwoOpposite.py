import h5py as h5
import numpy as np
import sys
from matplotlib import pyplot as plt

'''
@Brief: This script takes the simulated input from one waveguide, and then
    forms the mirror and then shift the two
'''
def main( argv ):
    d = 0.0
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print ("Usage: python crateTwoOpposite.py --d=<distance> --file=<farFieldH5file>")
            print ("d: Distance to shift the mirror in nm")
            print ("file: HDF5 file with the far field data")
        else:
            print ("Unknown argument %s"%(arg))
            return

    if ( fname == "" ):
        print ("No file specified!")
        return

    with h5.File( fname, 'r' ) as hf:
        dset = hf.get("exitIntensity")
        farField = hf.get("farField")
        attrs = dset.attrs
        wavenumber = float( farField.attrs["wavenumber"] )
        xmin = float( dset.attrs.get("xmin"))
        xmax = float( dset.attrs.get("xmax"))
        amp = np.array(dset)
        phase = np.array( hf.get("exitPhase") )
        #uid = int( attrs.get("uid") )
        farField = np.array( farField )
    field = amp*np.exp(1j*phase)
    mirror = field[-1::-1]
    dx = (xmax-xmin)/len(field)
    d=0.0
    happy = False
    while ( not happy ):
        dIndx = int(d/dx)
        newLength = len(field)+dIndx
        newField = np.zeros(newLength) + 1j*np.zeros(newLength)
        newField[:len(field)] = field
        newField[-len(field):] += mirror
        x = np.linspace(xmin,newLength*dx, newLength)
        plt.plot(x, np.abs(newField))
        plt.xlabel("x (nm)")
        plt.show( block=False )
        d = raw_input("Separation (nm) (type any non-numeric character to quit): ")
        try:
            d = float(d)
        except:
            happy = True
        plt.clf()


    fname = fname.split(".")[0]+"mirror.h5"
    farField = np.fft.fft( newField )
    np.fft.fftshift( farField )
    with h5.File(fname, 'w') as hf:
        dset = hf.create_dataset("exitIntensity", data=np.abs(newField))
        dset2 = hf.create_dataset("exitPhase", data=np.angle(newField))
        dset3 = hf.create_dataset("farField", data=np.array(farField) )
        dset.attrs["xmin"] = xmin
        dset.attrs["xmax"] = xmax
        dset.attrs["displacement"] = d
        dset3.attrs["wavenumber"] = wavenumber

    print ("New data written to %s"%(fname))

if __name__ == "__main__":
    main( sys.argv[1:])
