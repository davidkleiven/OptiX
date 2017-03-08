import sys
import json
import h5py as h5
import os

'''
Merges the HDF5 files created by one MPI run to one file
The name of the files is assumed to be
uid<UID>-dataname-otherstuff.h5

In addition it is assumed to be a JSON file describing the parameters with the name
uid<UID>-parameters.json
'''

def getDatasetName( fname ):
    datasetName = fname.split("-")[1]
    return datasetName

def getUID( fname ):
    return name.split("-")[0][3:]

def main( argv ):
    if ( len(argv) != 2 ):
        print ("Usage: python3 mergeMPIHDF5.py <uid> <folder>")
        return
    fname = argv[0]
    folder = argv[1]

    uid = getUID( fname )
    # Create a new hdf5 file
    outfname = "mergedMPI_"+uid+".h5"
    outfile = h5.File( outfname, 'w' )
    geom = outfile.create_group( "geometry" )
    flux = outfile.create_group( "flux" )

    # Parse the parameter file
    paramfname = "uid"+uid+"-parameters.json"
    paramFile = open( paramfname, 'r' )
    params = json.load( paramFile )
    paramFile.close()

    geom.attrs = params["geometry"]
    flux.attrs = params["dft"]

    # Loop through all files in the folder
    for filename in os.listdir(folder):
        if ( uid in filename ):
            name = getName( filename )
            if ( ".h5" in filename ):
                with h5.File( filename, 'r' ) as hf:
                    if ( "flux" in name ):
                        flux.create_dataset( name, np.array(hf.get(hf.keys()[0]) )
                    else:
                        geom.create_dataset( name, np.array(hf.get(hf.keys()[0])) )
            else if ( ".bin" in filename ):
                values = np.fromfile( filename, dtype=np.float64 )
                if ( "flux" in name ):
                    flux.create_dataset( name, values )
                else:
                    geom.create_dataset( name, values )
    outfile.close()
    print ( "Merged HDF5 file written to %s"%(outfname) )

if __name__ == "__main__":
    main( argv[1:] )
