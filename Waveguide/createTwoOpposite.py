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
        dx = farField.get("gridspacing")
        attrs = dset.attrs
        wavenumber = float( farField.attrs["wavenumber"] )
        ef = np.array( hf.get( "exitField" ) )
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
    dIndx = 0
    xmaxNew = xmax
    xminNew = xmin
    while ( not happy ):
        dIndx = int(d/dx)
        newLength = len(field)+np.abs(dIndx)
        center = int(newLength/2)
        newField = np.zeros(newLength) + 1j*np.zeros(newLength)
        if ( dIndx >= 0):
            newField[dIndx:len(field)+dIndx] = field
            newField[:len(field)] += mirror
            x = np.linspace(xmin,xmin+newLength*dx, newLength)
            xmaxNew = xmin+newLength*dx
            xminNew = xmin
        else:
            dIndx = np.abs(dIndx)
            newField[:len(field)] = field
            newField[dIndx:len(field)+dIndx] += mirror
            xmaxNew = xmax
            x = np.linspace(xmax-dx*newLength, xmax, newLength)
            xminNew = xmax-dx*newLength

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
    length = 32768
    padField = np.pad( newField, length/2, "edge")
    farField = np.fft.fft( padField )
    farField = np.fft.fftshift( farField )/np.sqrt(len(farField))
    with h5.File(fname, 'w') as hf:
        dset = hf.create_dataset("exitIntensity", data=np.abs(newField))
        dset2 = hf.create_dataset("exitPhase", data=np.angle(newField))
        dset3 = hf.create_dataset("farField", data=np.abs(farField) )
        dset4 = hf.create_dataset("exitField", data=np.real(newField))
        dsets = [dset, dset2, dset3, dset4]
        for ds in dsets:
            ds.attrs["xmin"] = xminNew
            ds.attrs["xmax"] = xmaxNew
            ds.attrs["displacement"] = d
            ds.attrs["wavenumber"] = wavenumber
            ds.attrs["gridspaceing"] = dx

    print ("New data written to %s"%(fname))

if __name__ == "__main__":
    main( sys.argv[1:])
