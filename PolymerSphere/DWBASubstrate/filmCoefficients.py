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
    print eps2
    mu1 = 1.0
    mu2 = 1.0
    incAngleTrans = np.arcsin( n1*np.sin(incAngle)/n2 )
    if ( pol == "TE" ):
        r = fe.rs( eps1, eps2, mu1, mu2, incAngle*180.0/np.pi, k ) 
        t01 = fe.ts( eps1, eps2, mu1, mu2, incAngle*180.0/np.pi, k)
        t10 = fe.ts(eps2, eps1, mu2, mu1, incAngleTrans*180.0/np.pi, k)
    else:
        r = fe.rp( eps1, eps2, mu1, mu2, incAngle*180.0/np.pi, k ) 
        t01 = fe.tp( eps1, eps2, mu1, mu2, incAngle*180.0/np.pi, k)
        t10 = fe.tp(eps2, eps1, mu2, mu1, incAngleTrans*180.0/np.pi, k)
    tfilm = (1.0-r**2)*np.exp(delta*1j)/( 1.0 - r**2 *np.exp(2.0*1j*delta))
    rfilm = r*(1.0-np.exp(2.0*1j*delta))/(1.0 - r**2*np.exp(2.0*1j*delta))
    return rfilm, tfilm
    
def main(argv):
    HELP_MSG = "Usage: python filmCoefficients-py --n2=<a+bj> --kd=10.0 [--help, --ksweep]\n"
    HELP_MSG += "help: Show this message\n"
    HELP_MSG += "n2: Refractive index of the film\n"
    HELP_MSG += "kd: Value of wave number times film thickness\n"
    HELP_MSG += "ksweep: Create plot from kd=0 to kd=10 for alpha = 0.5alpha_c, alpha_c and 2alpha_c\n"
    n1 = 1.0
    n2Found = False
    kdFound = False
    plotAngleSweep = True
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
        elif ( arg.find("--ksweep") != -1 ):
            plotAngleSweep = False
            print ("Running k-sweep instead....")
        else:
            print ("Unknown argument %s"%(arg))
            return 1

    if ( not n2Found ):
        print ("Did not find refrective index of film. Run python filmCoefficients.py --help")
        return 1
    elif ( (not kdFound) and plotAngleSweep ):
        print ("Did not find kd. Run with --help option\n")
        return 1

    d = 1.0
    if ( plotAngleSweep ):
        alphaMax = 5.0
        alpha_c = fe.criticalGrazingAngle(n1,n2)
            
        alpha = np.linspace(0.0, alphaMax, 10000.0)
        kd = 1E-2
        r1, t1 = filmCoefficients( n1, n2, d, kd, alpha )
        kd = 10
        r2, t2 = filmCoefficients( n1, n2, d, kd, alpha )
        kd = 1E4
        r3, t3 = filmCoefficients( n1, n2, d, kd, alpha )
        
        
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( alpha, np.abs(r1)**2, color=cs.COLORS[1], label="$kd=1$" )
        ax.plot( alpha, np.abs(t1)**2, color=cs.COLORS[1], ls="--" )
        ax.plot( alpha, np.abs(r2)**2, color=cs.COLORS[3], label="$kd=10^{2}$" )
        ax.plot( alpha, np.abs(t2)**2, color=cs.COLORS[3], ls="--" )
        ax.plot( alpha, np.abs(r3)**2, color=cs.COLORS[5], label="$kd=10^{4}$" )
        ax.plot( alpha, np.abs(t3)**2, color=cs.COLORS[5], ls="--" )
        ax.set_xlabel("$\\alpha/\\alpha_c$")
        ax.set_ylabel("Reflectance/Transmittance")
        ax.set_ylim(1E-6,9)
        ax.text(1.5, 1.5, "R (solid)")
        ax.text(2.5, 1.5, "T (dashed)")
        #ax.set_ylim(top=5)
        ax.set_yscale("log")
        ax.legend(loc="center right", frameon=False)
        fname = "Figures/filmCoeff.pdf"
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))
    else:
        k = np.logspace(-2, 5, 1000)
        incAngle = 0.5
        r1, t1 = filmCoefficients( n1, n2, d, k, incAngle )
        incAngle = 1.0
        r2, t2 = filmCoefficients( n1, n2, d, k, incAngle )
        incAngle = 2.0
        r3, t3 = filmCoefficients( n1, n2, d, k, incAngle )

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( k*d, np.abs(r1), color=cs.COLORS[1], label="$\\alpha = 0.5\\alpha_c$")
        ax.plot( k*d, np.abs(t1), color=cs.COLORS[1], ls="--")
        ax.plot( k*d, np.abs(r2), color=cs.COLORS[3], label="$\\alpha = \\alpha_c$")
        ax.plot( k*d, np.abs(t2), color=cs.COLORS[3], ls="--")
        ax.plot( k*d, np.abs(r3), color=cs.COLORS[5], label="$\\alpha = 2\\alpha_c$")
        ax.plot( k*d, np.abs(t3), color=cs.COLORS[5], ls="--")
        ax.set_xlabel("$kd$")
        ax.set_ylabel("Reflectance/Transmittance")
        ax.set_yscale("log")
        ax.set_ylim(1E-7,5.0)
        ax.set_xscale("log")
        ax.text(5E-2, 5E-2, "R (solid)\n T (dashed)")
        ax.legend(loc="lower center", frameon=False)
        fname = "Figures/fieldKSweep.pdf"
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))
         
if __name__ == "__main__":
    main(sys.argv[1:])
