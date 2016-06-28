import numpy as np
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json
import sys

FOLDERS_s = ["dataPlane/MultInc5s/WithEps", "dataPlane/MultInc20s/WithEps", "dataPlane/MultInc45s/WithEps", \
"dataPlane/MultInc75s/WithEps", "dataPlane/MultInc85s/WithEps"]
#FOLDERS_s = ["dataPlane/MultInc20/WithEps"]
FOLDERS_p = ["dataPlane/MultInc5p/WithEps", "dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", \
"dataPlane/MultInc75p/WithEps", "dataPlane/MultInc85p/WithEps"]
#FOLDERS_p = ["dataPlane/MultInc20p/WithEps"]
POLARISATRIONS = ["s", "p"]
SUBDIR="WithEps"
def transmittedTheta( theta, n1, n2 ):
    return np.arccos( np.sqrt( 1.0 - (n1*np.sin(theta)/n2)**2 ))

def rs(theta, n1, n2):
    theta = theta*np.pi/180.0
    thetaT = transmittedTheta(theta, n1, n2)
    return ( n1*np.cos(theta) - n2*np.cos(thetaT) )/( n1*np.cos(theta) + n2*np.cos(thetaT) )

def rp(theta, n1, n2): 
    theta = theta*np.pi/180.0
    thetaT = transmittedTheta(theta, n1, n2)
    return ( n1*np.cos(thetaT) - n2*np.cos(theta) )/( n1*np.cos(thetaT) + n2*np.cos(theta) )

def ts(theta, n1, n2):
    return 1.0+rs(theta, n1, n2)

def tp(theta, n1, n2):
    theta = theta*np.pi/180.0
    thetaT = transmittedTheta(theta, n1, n2)
    return 2.0*n1*np.cos(theta)/( n1*np.cos(thetaT) + n2*np.cos(theta) )
    
def main(argv):
    if ( len( argv ) < 2 ):
        print "[plotFieldCoeff] Usage: python plotFieldCoeff.py <ddirbase> <figure dir> <incident angles>"
        print "[plotFieldCoeff] Angle, polarisation and WithEps will be appended to the ddirbase"
        print "[plotFieldCoeff] Example: ddirbase = dataplane/MultInc --->"
        print "[plotFieldCoeff] 20 deg incidence s polarisation is located i dataplane/MultInc20s/WithEps"
        return
    ddir = argv[0]
    fdir = argv[1]
    try:
        incidentAngles = np.array(argv[2:]).astype(np.int32)
    except:
        print "[plotFieldCoeff] Error when converting incident angles to int..."
        return

    fig = plt.figure()
    ax = fig.add_subplot(111)
    figError = plt.figure()
    axError = figError.add_subplot(111)
    theta = np.linspace(0.0, 90.0, 101)
    n1 = 1.0
    n2 = 1.5
    ax.plot( theta, np.abs(rs(theta, n1, n2))**2, color='black' ) 
    ax.plot( theta, np.abs(rp(theta, n1, n2))**2, color='black' ) 
    ax.plot( theta, ts(theta, n1, n2)**2, color='black' )
    ax.plot( theta, tp(theta, n1, n2)**2, color='black' )

    msize = 2
    step = 10
    
    for pol in POLARISATRIONS:
        hasLabel = False
        for angle in incidentAngles:
            folder = ddir+"%d%s/WithEps"%(angle, pol)
            fname = folder+"/realFieldFourier.json"
            try:
                infile = open( fname, 'r')
                data = json.load(infile)
                infile.close()
            except:
                print ("Could not open file %s"%(fname))
                continue  

            reflectionNorm = np.array( data["reflected"]["norm"] )[0:-1:step]
            transmissionNorm = np.array( data["transmitted"]["norm"] )[0:-1:step]
            angleR = np.array( data["reflected"]["angle"] )[0:-1:step]
            angleT = np.array( data["transmitted"]["angle"] )[0:-1:step]
            errorRs = np.abs( reflectionNorm**2 - rs(angleR, n1, n2)**2 )/rs(angleR, n1, n2)**2
            errorTs = np.abs( transmissionNorm**2 - ts(angleT, n1, n2)**2)/ts(angleT, n1, n2)**2

            
            if ( pol == "s" ):
                markerR = 'o'
                markerT = 'x'
            else:
                markerR = 'x'
                markerT = '.'
            if ( not hasLabel ):
                ax.plot( angleR, reflectionNorm**2, markerR, color='black', ms=msize, fillstyle="none", \
                label="$|r_\mathrm{%s}|^2$"%(pol))
                ax.plot( angleT, transmissionNorm**2, markerT, color='black', ms=msize, label="$|t_\mathrm{%s}|^2$"%(pol))
                
                axError.plot( angleR, errorRs, markerR, color='black', ms=msize, fillstyle="none", label="$|r_\mathrm{%s}|^2$"%(pol) )
                axError.plot( angleT, errorTs, markerT, color='black', ms=msize, fillstyle="none", label="$|t_\mathrm{%s}|^2$"%(pol) )
                hasLabel = True
            else:
                ax.plot( angleR, reflectionNorm**2, markerR, color='black', ms=msize, fillstyle="none")
                ax.plot( angleT, transmissionNorm**2, markerT, color='black', ms=msize)

                axError.plot( angleR, errorRs, markerR, color='black', ms=msize, fillstyle="none")
                axError.plot( angleT, errorTs, markerT, color='black', ms=msize, fillstyle="none")

    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Incident angle (deg)")
    ax.set_ylabel("$|r_\mathrm{s}|^2$, $|r_\mathrm{p}|^2$, $|t_\mathrm{s}|^2$, $|t_\mathrm{p}|^2$")
    ax.legend(loc="upper left", ncol=2, frameon=False)
    axError.set_yscale('log')
    axError.set_xlabel("Incident angle (deg)")
    axError.set_ylabel("Relative error")
    axError.legend(loc='upper left', frameon=False)
    fnameCoeff = fdir+"/fieldCoefficients.pdf"
    fnameError = fdir+"/fieldCoefficientsError.pdf"
    fig.savefig( fnameCoeff, bbox_inches="tight" )
    figError.savefig( fnameError, bbox_inches="tight")
    print "[plotFieldCoeff] Figure written to %s"%(fnameCoeff)
    print "[plotFieldCoeff] Figure written to %s"%(fnameError)

if __name__ == '__main__':
    main(sys.argv[1:])
