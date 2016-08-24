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
import fresnelExact as fe 

def plotError( angle, error_rs, error_rp, error_ts, error_tp, brewster, fname ):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot( angle, error_rs, 'o', ms=2, color='black', fillstyle='none', label="$R_\mathrm{s}$" )
    ax.plot( angle, error_rp, 'x', ms=2, color='black', label="$R_\mathrm{p}$")
    if ( not error_ts is None ):
        ax.plot( angle, error_ts, '^', ms=2, color="black", label="$T_\mathrm{s}$")
    if ( not error_tp is None ):
        ax.plot( angle, error_tp, '.', ms=4, color="black", label="$T_\mathrm{p}$")
    ax.set_ylim( 1E-11, 1E-1 )
    if ( np.min(angle) < brewster ):
        ax.axvline( brewster, color='black', ls="dotted" )
        ax.text( brewster+1, 1E-8, "Brewster", rotation=-90 )
    ax.set_xlabel( "Incident angle (deg)" )
    ax.set_ylabel( "Relative error" )
    ax.set_yscale('log')
    ax.legend( loc='lower right', frameon=False )
    fig.savefig( fname, bbox_inches="tight" )
    print ( "Figure written to %s"%( fname ))

def main(argv):
    if ( len(argv) != 2 ):
        print ("Usage: python fluxPlot.py <datafile.json> --absorption=<true ro false>")
        return

    fig = plt.figure()
    figError = plt.figure()
    ax = fig.add_subplot(111)
    axError = figError.add_subplot(111)
    useAbsorption = False
    if ( argv[1] == "--absorption=true" ):
        useAbsorption = True
        print ("Using complex refractive indices...")
    else:
        print ("Using real refractive indices...")
    
    try:
        infile = open( argv[0], 'r' )
    except:
        print ("Error when opening file %s"%(argv[0]))
        return

    data = json.load(infile)
    infile.close()
    try:
        n2 = np.sqrt( float(data["geometry"]["EpsilonLow"]) )
        n1 = np.sqrt( float(data["geometry"]["EpsilonHigh"]) )
    except:
        n1 = 1.0
        n2 = 1.5

    incangle = np.array( data["IncidentAngle"] )
    angle = np.linspace(np.min(incangle)-0.05, 90.0, 1001)
    ax.plot( 90.0-np.array(data["IncidentAngle"]), data["FluxReflected"]["s"], 'o', ms=3, color='black', fillstyle='none', label="$R_\mathrm{TE}$")
    ax.plot( 90.0-np.array(data["IncidentAngle"]), data["FluxReflected"]["p"], 'x', ms=3, color='black', label="$R_\mathrm{TM}$" )
    if ( useAbsorption ):
        eps1 = 1.0
        eps2 = 0.99998*(1.0+2E-6*1j )
        mu1 = 1.0
        mu2 = 1.0
        k = 0.01
        ax.plot(90.0-angle, fe.Rs(eps1, eps2, mu1, mu2, angle, k), color='black')
        ax.plot(90.0-angle, fe.Rp(eps1, eps2, mu1, mu2, angle, k), color='black')
        ax.set_ylabel("Reflectance")
    else:
        ax.plot( 90.0-np.array(data["IncidentAngle"]), data["FluxTransmitted"]["s"], '^', ms=2, color='black', label="$T_\mathrm{s}$")
        ax.plot( 90.0-np.array(data["IncidentAngle"]), data["FluxTransmitted"]["p"], '.', ms=4, color='black', label="$T_\mathrm{p}$")
        ax.plot(90.0-angle, pf.Rs(90.0-angle, n1, n2), color='black')
        ax.plot(90.0-angle, pf.Rp(90.0-angle, n1, n2), color='black')
        ax.plot(90.0-angle, pf.Tp(90.0-angle, n1, n2), color='black')
        ax.plot(90.0-angle, pf.Ts(90.0-angle, n1, n2), color='black')
        ax.set_ylabel("Transmittance/Reflectance")
    ax.set_xlabel("Grazing incident angle (deg)")
    ax.legend( loc='center left', frameon=False)
    fname = "Figures/powerCoefficients.pdf"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

    # Plot relative error
    simulatedAngles = np.array( data["IncidentAngle"] )
    brewster = cra.brewster(n1,n2)*180.0/np.pi
    if ( useAbsorption ):
        Rs = fe.Rs(eps1,eps2,mu1,mu2,simulatedAngles,k)
        Rp = fe.Rp(eps1,eps2,mu1,mu2,simulatedAngles,k)
        Ts = None
        Tp = None
    else:
        Rs = pf.Rs(simulatedAngles, n1, n2)
        Rp = pf.Rp(simulatedAngles, n1, n2)
        Ts = pf.Ts(simulatedAngles, n1, n2)
        Tp = pf.Tp(simulatedAngles, n1, n2)
    errorRs = np.abs( np.array(data["FluxReflected"]["s"]) - Rs)/Rs
    errorRp = np.abs(np.array(data["FluxReflected"]["p"]) - Rp)/Rp
    if ( Ts is None ):
        errorTs = None
    else:
        errorTs = np.abs( np.array(data["FluxTransmitted"]["s"]) - Ts )/Ts
    if ( Tp is None ):
        errorTp = None
    else:
        errorTp = np.abs( np.array(data["FluxTransmitted"]["p"]) - Tp )/Tp
    plotError( simulatedAngles, errorRs, errorRp, errorTs, errorTp, brewster, "Figures/powerCoefficientsError.pdf" )
    
if __name__ == '__main__':
    main(sys.argv[1:] )
    
