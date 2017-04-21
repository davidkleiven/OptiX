import sys
import numpy as np
from scipy import optimize as opt
import json
import farFieldExactSphere as ffsph
import h5py as h5
from scipy import ndimage
from matplotlib import pyplot as plt

class CoatingConstrastOptimizer:
    def __init__( self ):
        self.coatedSphere = None
        self.qmin = 0.0
        self.qmax = 0.0
        self.ffSlice = None


    def costFunction( self, x ):
        assert( not self.coatedSphere is None )
        assert( not self.ffSlice is None )

        # Set the second refractive index to x
        self.coatedSphere.delta[1] = x
        q = np.linspace( self.qmin, self.qmax, len(self.ffSlice))
        formfactor = self.coatedSphere.formFactor(q)**2

        # Normalize
        formfactor *= np.sum(self.ffSlice)/np.sum(formfactor)
        return np.sum( (self.ffSlice - formfactor)**2 )

def main( argv ):
    if ( len(argv) != 1 ):
        print ( "Usage: python fitRefractiveindex.py inpufile.json" )

    try:
        optimizer = CoatingConstrastOptimizer()
        infile = open( argv[0] )
        params = json.load( infile )
        infile.close()

        optimizer.coatedSphere = ffsph.LayeredSphere()
        optimizer.coatedSphere.radii.append( params["R1"] )
        optimizer.coatedSphere.radii.append( params["R2"] )
        optimizer.coatedSphere.delta.append( 1.0 )
        optimizer.coatedSphere.delta.append( 1.0 )

        # Read data from h5 file
        with h5.File( params["datafile"], 'r' ) as hf:
            if ( not "/data/farField" in hf.keys() ):
                raise Exception("The hdf5 file needs to have a field named /data/farFeld")
            ff = hf.get("/data/farField")
            optimizer.qmin = ff.attrs.get("qmin")
            optimizer.qmax = ff.attrs.get("qmax")
            ffdata = np.array( ff )

        # Compute center of mass of the far field
        center = ndimage.center_of_mass( ffdata )
        optimizer.ffSlice = ffdata[center[0],:]
        del ffdata

        # Optimize the ratio
        optimalRatio = opt.minimize( optimizer.costFunction, 5.5 )

        # Set the density to the optimal value
        optimizer.coatedSphere.delta[1] = optimalRatio.x

        print ("Optimal density ratio: %.4E"%(optimalRatio.x))

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        q = np.linspace( optimizer.qmin, optimizer.qmax, len(optimizer.ffSlice) )
        ffactor = optimizer.coatedSphere.formFactor(q)**2
        ffactor *= np.max(optimizer.ffSlice)/np.max(ffactor)
        ax.plot(  q, optimizer.ffSlice, color="blue")
        ax.plot( q, ffactor, color="red")
        ax.set_yscale("log")
        plt.show()
    except Exception as exc:
        print (str(exc))
        return 1
    return 0

if __name__ == "__main__":
    main( sys.argv[1:] )
