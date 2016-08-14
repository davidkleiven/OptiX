import sys
sys.path.append( "../FresnelFDTD/" )
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update( mplLaTeX.params )
from matplotlib import pyplot as plt
import numpy as np

OMEGA = 2.0*np.pi*0.5
def rhs( ky, epsscat ):
    return np.sqrt( OMEGA**2 *(1.0-epsscat)/ky**2 - 1.0 )
def main(argv):
    if ( len(argv) != 2 ):
        print ("Usage: python plotTransientEq.py --fig=<figname> --epscladding=<eps in cladding>")
        return 1

    # Parse arguments
    figname = ""
    epsclad = -1.0
    for arg in argv:
        if ( arg.find("--fig=") != -1 ):
           figname = arg.split("--fig=")[1]
        elif ( arg.find("--epscladding=") != -1 ):
            epsclad = float( arg.split("--epscladding=")[1] )
    
    # Consistency check
    if ( figname == "" ):
        print ("No figure name specified...")
        return 1
    elif ( epsclad < 0.0 ):
        print ("No epsilon in cladding specified...") 
        return 1
    ky = np.linspace(0.0, 2.0*np.pi/2.0, 200 )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot( ky, np.tan(ky), 'k', label="$\\tan{k_y}$")
    ax.plot( ky, rhs(ky, epsclad), 'k--', label="$\sqrt{\\frac{\omega^2 (\epsilon_2 - \epsilon_1)}{k_y^2} - 1 }$")
    ax.set_ylim(0.0, 4.0)
    ax.set_xlabel( "$k_y$" )
    ax.legend(loc="upper right", frameon=False)
    fig.savefig(figname, bbox_inches="tight" )
    print ("Figure written to %s"%(figname))

if __name__ == "__main__":
    main(sys.argv[1:])
    
