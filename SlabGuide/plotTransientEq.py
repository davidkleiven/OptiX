import sys
sys.path.append( "../FresnelFDTD/" )
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update( mplLaTeX.params )
from matplotlib import pyplot as plt
import numpy as np
from scipy import optimize as opt

OMEGA = 2.0*np.pi*1.0
def rhs( ky, epsscat ):
    return np.sqrt( OMEGA**2 *(1.0-epsscat)/ky**2 - 1.0 )

def solverFunc( ky, epsscat ):
    return np.tan(ky) - rhs(ky,epsscat)

def main(argv):
    global OMEGA
    if ( len(argv) != 4 ):
        print ("Usage: python plotTransientEq.py --fig=<figname> --epscladding=<eps in cladding> --freq=<frequency>")
        print ("--guess=<unity or lowfreq>\n")
        return 1

    # Parse arguments
    figname = ""
    epsclad = -1.0
    initguess = ""
    for arg in argv:
        if ( arg.find("--fig=") != -1 ):
           figname = arg.split("--fig=")[1]
        elif ( arg.find("--epscladding=") != -1 ):
            epsclad = float( arg.split("--epscladding=")[1] )
        elif ( arg.find("--freq=") != -1 ):
            OMEGA *= float( arg.split("--freq=")[1] )
        elif ( arg.find("--guess=") ):
            initguess = arg.split("--guess=")[1]
    
    # Consistency check
    if ( figname == "" ):
        print ("No figure name specified...")
        return 1
    elif ( epsclad < 0.0 ):
        print ("No epsilon in cladding specified...") 
        return 1
    elif ( initguess == "" ):
        print ("Did not find any init gues. Using lowfrew...")
        initguess = "lowfreq"
    kymax = OMEGA*np.sqrt(1.0-epsclad)
    if ( kymax > np.pi/2.0):
        kymax = np.pi/2.0
    ky = np.linspace(1E-10, kymax, 2000 )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot( ky, np.tan(ky), 'k', label="$\\tan{k_y}$")
    ax.plot( ky, rhs(ky, epsclad), 'k--', label="$\sqrt{\\frac{\omega^2 (\epsilon_2 - \epsilon_1)}{k_y^2} - 1 }$")
    ax.set_ylim(0.0, 4.0)
    ax.set_xlabel( "$k_y$" )
    ax.legend(loc="upper right", frameon=False)
    fig.savefig(figname, bbox_inches="tight" )
    print ("Figure written to %s"%(figname))

    # Solve the equation
    if ( initguess == "lowfreq" ):
        guess = OMEGA**2 *(1.0-epsclad)
    else:
        guess = 1.0
    try:
        kyAnswer = opt.newton( solverFunc, guess, args=(epsclad,), maxiter=1000)
    except:
        print ("Failed to converge, using small ky expansion")
        kyAnswer = OMEGA*np.sqrt(1.0-epsclad)
    print ("Intersection point %.4E"%(kyAnswer))
    

if __name__ == "__main__":
    main(sys.argv[1:])
    
