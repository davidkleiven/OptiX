import sys
sys.path.append("Scripts")
import numpy as np
from matplotlib import pyplot as plt
import harminvSingle as hinv
import h5py as h5

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python bandDiagram.py dataset.h5")
        return 1

    plotter = hinv.BlochRun()
    with h5.File(argv[0],'r') as hf:
        plotter.buildBlochRegs( hf )
    plotter.plotBands( maxbands=6 )
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
