import sys
sys.path.append("../../FresnelFDTD")
import numpy as np
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
from matplotlib import pyplot as plt

def main():
    delta = 4.14E-5
    beta = 3.45E-6

    n2 = 1 - delta + 1j*beta
    ac = np.sqrt(2.0*delta)
    angle = np.linspace(0.0, 2.0, 101)*ac
    angle_t = np.real( np.sqrt(1.0 - (np.cos(angle)/n2)**2) )

    print np.real( np.sqrt(( 1- np.cos(0.0)**2 * n2**(-2))) )
    asymptote = beta/ac**2 *(1 + 0.5*(angle/ac)**2)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot ( angle/ac, angle_t/ac, color="black")
    ax.plot( angle/ac, asymptote, color="black", ls="--")
    #ax.axhline(asymptote, ls="--",color="black")
    #ax.text( 0.8*angle[-1]/ac, asymptote, "$\\frac{\\beta}{\\alpha_c}$")
    ax.set_xlabel("$\\frac{\\alpha}{\\alpha_c}$")
    ax.set_ylabel("$\\frac{\\alpha_t}{\\alpha_c}$", rotation=0)
    ax.set_yscale("log")
    fname = "../Figures/transAngle.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main()
