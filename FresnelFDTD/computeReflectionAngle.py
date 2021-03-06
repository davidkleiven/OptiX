#coding: utf-8
import numpy as np
import mplLaTeX
import matplotlib
matplotlib.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
from scipy import stats
from scipy import optimize
import json
import sys
ANGLES = [5, 20, 45, 75, 85]
POLARISATIONS=["s","p"]
BASE = "dataPlane/MultInc"
SUBDIR = "WithEps"
FNAME = "realFieldFourier.json"

def findReflectionAngle(theta_r, theta_in, waveNumber, distanceFromSlab, phase):
    # Force the solver to stay within +- pi
    if (( theta_r > np.pi/2.0 ) or ( theta_r < 0.0 )):
        return np.inf
    return waveNumber*distanceFromSlab*(np.sin(theta_in)*(np.tan(theta_in)+np.tan(theta_r)) \
    - 1.0/np.cos(theta_in) - 1.0/np.cos(theta_r)) - phase

def expectedPhase(freq, fcenter, thetaCenter, distanceFromPlane):
    sqrtArg = freq**2 - (fcenter*np.sin(thetaCenter))**2 
    sqrtArg[sqrtArg < 0.0] = 0.0
    return 2.0*np.pi*2.0*distanceFromPlane*np.sqrt( sqrtArg )
 
def transmissionAngle( n1, n2, theta_in ):
    sinTheta = n1*np.sin(theta_in)/n2 
    sinTheta[ np.abs(sinTheta) > 1.0 ] = 1.0
    return np.arcsin( sinTheta)

def expectedTransmissionPathTime(n1, n2, theta_i, distanceYFromSlab):
    theta_t = transmissionAngle( n1, n2, theta_i )
    d = distanceYFromSlab*np.sqrt( np.cos(theta_t)**(-2) + np.cos(theta_i)**(-2) - \
                    2.0*np.cos(theta_i-theta_t)/(np.cos(theta_i)*np.cos(theta_t)) )
    return n1*d*np.sin(theta_i) + distanceYFromSlab*( n2/np.cos(theta_t) - n1/np.cos(theta_i) )

def brewster(n1, n2):
    return np.arctan(n2/n1)

def main(argv):
    if ( len( argv ) < 2 ):  
        print "[phase] Usage: python computeReflectionAngle.py <ddirbase> <figure dir> <incident angles>"
        print "[phase] Angle, polarisation and WithEps will be appended to the ddirbase"
        print "[phase] Example: ddirbase = dataplane/MultInc --->"
        print "[phase] 20 deg incidence s polarisation is located i dataplane/MultInc20s/WithEps"
        return
    ddir = argv[0]
    fdir = argv[1]
    try:
        incidentAngles = np.array(argv[2:]).astype(np.int32)
    except:
        print "[phase] Error when converting incident angles to int..."
        return
    fig = plt.figure()
    figSign = plt.figure()
    axSign = figSign.add_subplot(111)
    ax = fig.add_subplot(111)
    figSignT = plt.figure()
    axSignT = figSignT.add_subplot(111)

    figT = plt.figure()
    axT = figT.add_subplot(111)

    tanReflTimesTanAngle = None
    tanAngle = None
    step = 5
    brewsterAngle = -1.0
    distanceFromPlane = -1.0
    for pol in POLARISATIONS:
        hasLabel = False
        if ( pol == "s" ):
            marker = "o"
        else:
            marker = "x"
        for theta in ANGLES:
            fname = ddir+"%d%s/"%(theta, pol)+SUBDIR +"/"+FNAME
            try:
                infile = open(fname,'r')
                data = json.load(infile)
                infile.close()
            except:
                print ("[phase] Could not find file %s"%(fname))
                continue

            if ( brewsterAngle < 0.0 ):
                n1 = np.sqrt( data["geometry"]["EpsilonLow"])
                n2 = np.sqrt( data["geometry"]["EpsilonHigh"] )
                brewsterAngle = brewster(n1, n2)
            distanceFromPlane = data["reflected"]["position"] - data["geometry"]["slabPosition"]
            transmissionMonitorDistance = np.abs(data["transmitted"]["position"] - data["geometry"]["slabPosition"])
            phase = np.array( data["reflected"]["phase"] )
            angle = np.array( data["reflected"]["angle"] )*np.pi/180.0
            phaseTransmitted = np.array( data["transmitted"]["phase"] )
            angleTransmitted = np.array( data["transmitted"]["angle"] )*np.pi/180.0
            k = 2.0*np.pi*np.array( data["reflected"]["frequency"] )
            omegaTransmitted = 2.0*np.pi*np.array( data["transmitted"]["frequency"] )

            angle = angle[::step]
            k = k[::step]
            phase = phase[::step]
            angleTransmitted = angleTransmitted[::step]
            omegaTransmitted = omegaTransmitted[::step]
            phaseTransmitted = phaseTransmitted[::step]
            kTransmitted = omegaTransmitted

            # Compute expected phase for the reflected fields
            x = 2.0*distanceFromPlane*k*np.cos(angle)
            mNoSignChange =  (-x-phase)/(2.0*np.pi)
            mSignChange = (-x+np.pi-phase)/(2.0*np.pi)
            m = np.zeros(len(mSignChange))
            for i in range(0, len(mNoSignChange)):
                diffSign = np.abs( np.round(mSignChange[i]) - mSignChange[i]) 
                diffNoSign = np.abs( np.round(mNoSignChange[i]) - mNoSignChange[i] )
                if ( diffNoSign < diffSign ):
                    m[i] = np.round(mNoSignChange[i])
                else:
                    m[i] = np.round(mSignChange[i])
            phase += m*2.0*np.pi

            # Compute expected phase for the transmitted fields
            theta_t = transmissionAngle(n1,n2,angleTransmitted) 
            xTransmitted = omegaTransmitted*expectedTransmissionPathTime( n1, n2, angleTransmitted, transmissionMonitorDistance )
            label="$\\angle r_\mathrm{s}$"
            tlabel="$\\angle t_\mathrm{s}$"
            msize = 5
        
            m = np.round( (-xTransmitted - phaseTransmitted)/(2.0*np.pi) )
            phaseTransmitted += 2.0*np.pi*m
            # Sign of the reflected
            signPhase = phase + 2.0*k*distanceFromPlane*np.cos(angle)
            signPhaseT = phaseTransmitted + xTransmitted
            if ( pol == "p" ):
                msize=2
                label="$\\angle r_\mathrm{p}$"
                tlabel="$\\angle t_\mathrm{p}$"
            if ( hasLabel ):
                ax.plot(x, phase, marker, color="black", markersize=msize, fillstyle="none")
                axSign.plot( angle*180.0/np.pi, signPhase, marker, color='black', markersize=msize, fillstyle="none")
                axT.plot( xTransmitted, phaseTransmitted, marker, color='black', markersize=msize, fillstyle="none")
                axSignT.plot( angleTransmitted*180.0/np.pi, signPhaseT, marker, color='black', markersize=msize, fillstyle="none")
            else: 
                ax.plot(x, phase, marker, color="black", markersize=msize, fillstyle="none", label=label)
                axSign.plot( angle*180.0/np.pi, signPhase, marker, color='black', markersize=msize, fillstyle="none", label=label)
                axT.plot( xTransmitted, phaseTransmitted, marker, color='black', markersize=msize, fillstyle="none", label=tlabel)
                axSignT.plot( angleTransmitted*180.0/np.pi, signPhaseT, marker, color='black', markersize=msize, \
                                fillstyle="none", label=tlabel)
                hasLabel = True

    x1, x2 = ax.get_xlim()
    x = np.linspace(0.9*x1, 1.1*x2, 11)
    ax.plot( x, -x+np.pi, '--', color='black', label="Sign change")
    ax.plot( x, -x, color='black', label="No sign change")
    ax.set_xlabel("$\phi_\mathrm{path}=\\frac{2y\omega}{c} \cos \\theta_\mathrm{i}$")
    ax.set_ylabel("$\\phi_\mathrm{\omega}$", rotation=0)

    # Set labels in multiples of pi
    y1, y2 = ax.get_ylim()
    mMax = int(np.floor( y2/np.pi))
    mMin = int(np.floor( y1/np.pi))
    yticks = np.arange(mMin, mMax+1)*np.pi
    x1, x2 = ax.get_xlim()
    mxMax = int(np.ceil(x2/np.pi))
    xticks = np.arange(0, mxMax)*np.pi
    ax.set_yticks(yticks)
    ax.set_xticks(xticks)
    ylabels = []
    for i in range(mMin, mMax+1):
        if ( i==0 ):
            ylabels.append("$0$")
        elif ( i==1 ):
            ylabels.append("$\pi$")
        elif ( i== -1):
            ylabels.append("$-\pi$")
        else:
            ylabels.append("$%d\pi$"%(i))
    ax.set_yticklabels(ylabels)
    xlabels = []
    for i in range(0, mxMax):
        if ( i==0):
            xlabels.append("$0$")
        elif (i==1):
            xlabels.append("$\pi$")
        else:
            xlabels.append("$%d\pi$"%(i))
    ax.set_xticklabels(xlabels)
    ax.legend(loc="upper right", frameon=False)
    fnameRefPhase = fdir + "/reflectedPhase.pdf"
    fig.savefig(fnameRefPhase, bbox_inches="tight")
    print "[phase] Figure written to %s"%(fnameRefPhase)

    # Fix the sign plot
    axSign.set_xlabel("Incident angle (deg)")
    axSign.set_ylabel("$\phi_\omega - \phi_\mathrm{path}$")
    axSign.legend(loc='lower left', frameon=False)
    axSign.set_yticks([-np.pi, -np.pi/2.0, 0.0, np.pi/2.0, np.pi])
    axSign.set_yticklabels(["$-\pi$", "$-\\frac{\pi}{2}$", "$0$", "$\\frac{\pi}{2}$", "$\pi$"])
    axSign.axvline(brewsterAngle*180.0/np.pi, color='black', ls='--')
    y1, y2 = axSign.get_ylim()
    axSign.text( brewsterAngle*180.0/np.pi, 0.3*(y2-y1)+y1, "Brewster", rotation=90)
    fname = fdir + "/reflectedSign.pdf"
    figSign.savefig(fname, bbox_inches="tight")
    print "[phase] Figure written to %s"%(fname)

    # Fix transmission plot
    axT.set_xlabel("$\\phi_{t,\mathrm{path}}$")
    axT.set_ylabel("$\\phi_{\mathrm{t},\omega}$")
    x1, x2 = axT.get_xlim()
    x = np.linspace(0.0, 1.1*x2, 11)
    axT.plot(x, -x, color='black')
    x1, x2 = axT.get_xlim()
    m = int(np.floor(x2/np.pi))
    xticks = np.arange(0, m+1)*np.pi
    axT.set_xticks(xticks)
    xlabels = []
    xlabels.append("$0$")
    xlabels.append("$\pi$")
    for i in range(2,m+1):
        xlabels.append("$%d\pi$"%(i))
    axT.set_xticklabels(xlabels)

    y1, y2 = axT.get_ylim()
    mMin = int(np.ceil(y1/np.pi))
    mMax = int(np.floor(y2/np.pi))
    yticks = np.arange(mMin, mMax)*np.pi
    ylabels = []
    for m in range(mMin, mMax+1):
        if ( m==0 ):
            ylabels.append("$0$")
        elif ( m == 1):
            ylabels.append("$\pi$")
        elif (m == -1 ):
            ylabels.append("$-\pi$")
        else:
            ylabels.append("$%d\pi$"%(m))
    axT.set_yticks(yticks)
    axT.set_yticklabels(ylabels)
    axT.legend(loc='upper right', frameon=False)
    fname = fdir + "/transmittedPhase.pdf"
    figT.savefig(fname, bbox_inches="tight")
    print "[phase] Figure written to %s"%(fname)

    # Transmission sign plot
    axSignT.set_xlabel("Incident angle (deg)")
    axSignT.set_ylabel("$\phi_{\mathrm{t},\omega} - \phi_{\mathrm{t},\mathrm{path}}$")
    axSignT.legend(loc='upper right', frameon=False)    
    axSignT.set_yticks([-np.pi/4.0, -np.pi/8.0, 0.0, np.pi/8.0, np.pi/4.0])
    axSignT.set_yticklabels(["$-\\frac{\pi}{4}$", "$-\\frac{\pi}{8}$", "$0$", "$\\frac{\pi}{8}$", "$\\frac{\pi}{4}$"])
    fname = fdir + "/transmissionSign.pdf"
    figSignT.savefig(fname, bbox_inches="tight")
    print "[phase] Figure written to %s"%(fname)
    
if __name__ == "__main__":
    main(sys.argv[1:])
