import sys
sys.path.append( "../FresnelFDTD" )
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
import numpy as np
from matplotlib import pyplot as plt
import json
import plotFieldCoeff as pfc

def main(argv):
    if ( len(argv) != 1 ):
        print ("Usage: python amplitudes.py <datafile.json>")
        return

    fig = plt.figure()
    ax = fig.add_subplot(111)

    try:
        infile = open( argv[0], 'r' )
    except:
        print ("Error when opening file %s"%(argv[0]))
        return
    data = json.load(infile)
    infile.close()

    try:
        n1 = np.sqrt( float(data["geometry"]["EpsilonHigh"]) )
        n2 = np.sqrt( float(data["geometry"]["EpsilonLow"]) )
    except:
        print ("Could not find EpsilonHigh and EpsilonLow in data file. Using defaults...")
        n1 = 1.0
        n2 = 1.5

    ax.plot( data["IncidentAngle"], data["AmplitudeReflected"]["s"], 'o', ms=2, color='black', fillstyle='none', \
    label="$|r_\mathrm{s}|^2$" )

    ax.plot( data["IncidentAngle"], data["AmplitudeReflected"]["p"], 'x', ms=2, color='black', label="$|r_\mathrm{p}|^2$" )
    ax.plot( data["IncidentAngle"], data["AmplitudeTransmitted"]["s"], '^', ms=2, color='black', label="$|t_\mathrm{s}|^2$" )
    #ax.plot( data["IncidentAngle"], data["AmplitudeTransmitted"]["p"], '.', ms=2, color='black', label="$|t_\mathrm{p}|^2$" )

    # Plot exact values
    angle = np.linspace(0.0, 90.0, 101)
    ax.plot( angle, pfc.rs(angle, n1, n2)**2, color='black')
    ax.plot( angle, pfc.rp(angle, n1, n2)**2, color='black')
    ax.plot( angle, pfc.ts(angle, n1, n2)**2, color='black')
    ax.plot( angle, pfc.tp(angle, n1, n2)**2, color='black')
    
    ax.set_xlabel("Incident angle")
    ax.set_ylabel("Reflection/Transmission amplitudes")
    ax.legend(loc='center left', frameon=False)
    fname = "Figures/amplitudes.pdf"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname ))
    
if __name__ == '__main__':
    main(sys.argv[1:])
