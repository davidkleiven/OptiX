import sys
sys.path.append("Scripts")
import plotSpectrum as pspec
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"]=False
from matplotlib import pyplot as plt

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python crossSectionWithChlorophyll.py <hdf5 file>")
        return

    spec = pspec.initSpectrum( argv[0] )
    try:
        spec.SrDistribution()
    except Exception as exc:
        print (str(exc))
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
