import sys
sys.path.append("../../FresnelBEM")
import fresnelExact as fe
sys.path.append("../../FresnelFDTD")
sys.path.append("../../")
import matplotlib
import mplLaTeX as ml
matplotlib.rcParams.update( ml.params )
import colorScheme as cs
import numpy as np
from matplotlib import pyplot as plt

def phaseDifference( n1, n2, thickness, k, incAngleRadians ):
    ''' 
    @brief: Computes the phase difference delta in the film coefficients.
            It enters through exp(2i*delta)
    '''
    kz_t_unit = np.sqrt( 1.0 - (n1*np.sin(incAngleRadians)/n2)**2 )
    return n2*k*thickness*kz_t_unit

def filmCoefficients( n1, n2, filmThick, k, incAngleCriticalAngle, pol="TE" ):
    '''
    @brief Computes the film reflection and transission coefficient
    @return reflection, transmission
    '''
    alpha_c = fe.criticalGrazingAngle( n1, n2 )
    alpha = alpha_c*incAngleCriticalAngle
    print "Critical angle: %.3f"%(alpha_c*180.0/np.pi)
    
    incAngle = np.pi/2.0-alpha
    
    delta = phaseDifference( n1, n2, filmThick, k, incAngle )
    eps1 = np.sqrt(n1)
    eps2 = np.sqrt(n2)
    mu1 = 1.0
    mu2 = 1.0
    if ( pol == "TE" ):
        r = fe.rs( eps1, eps2, mu1, mu2, incAngle*180.0/np.pi, k ) 
        t01 = fe.ts( eps1, eps2, mu1, mu2, incAngle*180.0/np.pi, k)
        t10 = fe.ts(eps2, eps1, mu2, mu1, incAngle*180.0/np.pi, k)
    else:
        r = fe.rp( eps1, eps2, mu1, mu2, incAngle*180.0/np.pi, k ) 
        t01 = fe.tp( eps1, eps2, mu1, mu2, incAngle*180.0/np.pi, k)
        t10 = fe.tp(eps2, eps1, mu2, mu1, incAngle*180.0/np.pi, k)
    tfilm = t01*t10*np.exp(delta*1j)/( 1.0 - r**2 *np.exp(2.0*1j*delta))
    return r*(1.0-np.exp(2.0*1j*delta)/(1.0-r**2 *np.exp(2.0*1j*delta))), tfilm
    
def main(argv):
    HELP_MSG = "Usage: python filmCoefficients-py --n2=<a+bj> --kd=<wave>[--help]\n"
    HELP_MSG += "help: Show this message\n"
    HELP_MSG += "n2: Refractive index of the film\n"
    HELP_MSG += "kd: Wave number times film thickness\n"
    n1 = 1.0
    n2Found = False
    kdFound = False
    for arg in argv:
        if ( arg.find("--n2=") != -1 ):
            n2 = complex(arg.split("--n2=")[1])
            n2Found = True
        elif ( arg.find("--help") != -1 ):
            print HELP_MSG
            return 0
        elif ( arg.find("--kd=") != -1 ):
            kd = float(arg.split("--kd=")[1])
            kdFound = True
        else:
            print ("Unknown argument %s"%(arg))
            return 1

    if ( not n2Found ):
        print ("Did not find refrective index of film. Run python filmCoefficients.py --help")
        return 1
    elif ( not kdFound ):
        print ("Did not find kd. Run with --help option\n")
        return 1

    d = 1.0
    alphaMax = 3.0
    alpha_c = fe.criticalGrazingAngle(n1,n2)
        
    alpha = np.linspace(0.0, alphaMax, 1000.0)
    r, t = filmCoefficients( n1, n2, d, kd, alpha )
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( alpha, np.abs(r), "k", label="$|r|$" )
    ax.plot( alpha, np.abs(t), "k--", label="$|t|$" )
    ax.set_xlabel("$\\alpha/\\alpha_c$")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(1E-5,1E2)
    ax.set_yscale("log")
    ax.legend(loc="upper right", frameon=False)
    fname = "Figures/filmCoeff.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))
if __name__ == "__main__":
    main(sys.argv[1:])
