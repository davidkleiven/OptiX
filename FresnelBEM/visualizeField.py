import sys
import numpy as np
from mayavi import mlab
import json

def main(argv):
    if ( len(argv) != 1 ):
        print ("Usage: python visualizeField.py <datafile.json>")
        return

    # Read input file
    try:
        infile = open( argv[0], 'r' )
    except:
        print ("Error when opening file %s"%(argv[0]))
        return
    data = json.load( infile )
    infile.close()

    y = np.array( data["points"]["y"] )
    z = np.array( data["points"]["z"] )
    x = np.zeros( len(y) )
    Ex = np.array( data["field"]["x"] )
    pts = mlab.points3d( y, z, x, Ex, scale_mode="scalar", scale_factor=0.0, mode="point")
    mesh = mlab.pipeline.delaunay2d( pts )
    surf = mlab.pipeline.surface( mesh )
    
    fname = "Figures/fieldyzPlane.png"
    mlab.view(0.0,0.0,1.0, (0.0,0.0,0.0))
    mlab.scalarbar()
    mlab.xlabel( "y" )
    mlab.ylabel( "z" )
    mlab.show()
    #mlab.savefig( fname )
    
if __name__ == '__main__':
    main( sys.argv[1:] )
