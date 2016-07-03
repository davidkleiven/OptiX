import sys
sys.path.append( "../FresnelFDTD/" )
import mplLaTeX
import plotFlux as pf # Contains the exact expression for the fluxes
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
import numpy as np
from matplotlib import pyplot as plt
import json
import computeReflectionAngle as cra

def main(argv):
    if ( len(argv) != 1 ):
        print ("Usage: python fluxPlot.py <datafile.json>")
        return

    fig = plt.figure()
    figError = plt.figure()
    ax = fig.add_subplot(111)
    axError = figError.add_subplot(111)
    
    try:
        infile = open( argv[0], 'r' )
    except:
        print ("Error when opening file %s"%(argv[0]))
        return

    data = json.load(infile)
    infile.close()
    try:
        n1 = np.sqrt( float(data["geometry"]["EpsilonLow"]) )
        n2 = np.sqrt( float(data["geometry"]["EpsilonHigh"]) )
    except:
        n1 = 1.0
        n2 = 1.5

    angle = np.linspace(0.0, 90.0, 101)
    ax.plot( data["IncidentAngle"], data["FluxReflected"]["s"], 'o', ms=2, color='black', fillstyle='none', label="$R_\mathrm{s}$")
    ax.plot( data["IncidentAngle"], data["FluxReflected"]["p"], 'x', ms=2, color='black', label="$R_\mathrm{p}$" )
    ax.plot( data["IncidentAngle"], data["FluxTransmitted"]["s"], '^', ms=2, color='black', label="$T_\mathrm{s}$")
    ax.plot( data["IncidentAngle"], data["FluxTransmitted"]["p"], '.', ms=4, color='black', label="$T_\mathrm{p}$")
    ax.plot(angle, pf.Rs(angle, n1, n2), color='black')
    ax.plot(angle, pf.Rp(angle, n1, n2), color='black')
    ax.plot(angle, pf.Tp(angle, n1, n2), color='black')
    ax.plot(angle, pf.Ts(angle, n1, n2), color='black')
    ax.set_xlabel("Incident angle")
    ax.set_ylabel("Transmittance/Reflectance")
    ax.legend( loc='center left', frameon=False)
    fname = "Figures/powerCoefficients.pdf"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

    # Plot relative error
    fig = plt.figure()
    ax = fig.add_subplot(111)
    simulatedAngles = np.array( data["IncidentAngle"] )
    brewster = cra.brewster(n1,n2)*180.0/np.pi
    Rs = pf.Rs(simulatedAngles, n1, n2)
    Rp = pf.Rp(simulatedAngles, n1, n2)
    Ts = pf.Ts(simulatedAngles, n1, n2)
    Tp = pf.Tp(simulatedAngles, n1, n2)
    errorRs = np.abs( np.array(data["FluxReflected"]["s"]) - Rs)/Rs
    errorRp = np.abs(np.array(data["FluxReflected"]["p"]) - Rp)/Rp
    errorTs = np.abs( np.array(data["FluxTransmitted"]["s"]) - Ts )/Ts
    errorTp = np.abs( np.array(data["FluxTransmitted"]["p"]) - Tp )/Tp
    ax.plot( simulatedAngles, errorRs, 'o', ms=2, color='black', fillstyle='none', label="$R_\mathrm{s}$" )
    ax.plot( simulatedAngles, errorRp, 'x', ms=2, color='black', label="$R_\mathrm{p}$")
    ax.plot( simulatedAngles, errorTs, '^', ms=2, color="black", label="$T_\mathrm{s}$")
    ax.plot( simulatedAngles, errorTp, '.', ms=4, color="black", label="$T_\mathrm{p}$")
    ax.set_ylim( 1E-11, 1E-1 )
    ax.axvline( brewster, color='black', ls="dotted" )
    ax.text( brewster+1, 1E-8, "Brewster", rotation=-90 )
    ax.set_xlabel( "Incident angle (deg)" )
    ax.set_ylabel( "Relative error" )
    ax.set_yscale('log')
    ax.legend( loc='lower right', frameon=False )
    fname = "Figures/powerCoefficientsError.pdf"
    fig.savefig( fname, bbox_inches="tight" )
    print ( "Figure written to %s"%( fname ))
if __name__ == '__main__':
    main(sys.argv[1:] )
    
