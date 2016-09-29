import grazing as graz
import grazingTransmission as gzt
import numpy as np
import matplotlib
import mplLaTeX as ml
matplotlib.rcParams.update( ml.params )
from matplotlib import pyplot as plt
import colorScheme as cs

class CombinedDetector:
    def __init__(self):
        self.n2 = 0.992+0.002j
        self.thick = 2.0
        self.reflected = graz.GrazingHandler(True)
        self.transmitted = gzt.GrazingTransmissionHandler()
        self.reflected.setEpsilonSubst(self.n2**2)
        self.transmitted.setEpsilonSubst(self.n2**2)
        self.reflected.setFilmThickness(self.thick)
        self.transmitted.setFilmThickness(self.thick)

    def plotTotal(self, angles):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for i in range(0, len(angles)):
            if ( angles[i] < 1.0 ):
                label = "$\\alpha = %.1f\\alpha_c$"%(angles[i])
            elif ( angles[i] == 1 ):
                label = "$\\alpha = \\alpha_c$"%(angles[i])
            else:
                label = "$\\alpha=%d\\alpha_c$"%(angles[i])
            self.reflected.prepareDWBA(angles[i])
            self.transmitted.prepareDWBA(angles[i]) 
            alpha_f_plus = np.arctan(self.reflected.x/self.reflected.detectorPosition)/self.reflected.alpha_c
            alpha_f_minus = np.arctan(self.transmitted.x/self.transmitted.detectorPosition)/self.transmitted.alpha_c

            
            totPluss = np.abs(self.reflected.bornTotal())**2
            firstBornPluss = np.abs(self.reflected.f1)**2
            firstBornMinus = np.abs(self.transmitted.born())**2
            totMinus = np.abs(self.transmitted.total())**2
        
            ax.plot( alpha_f_plus, totPluss, color=cs.COLORS[i], label=label)
            ax.plot( alpha_f_minus, totMinus, color=cs.COLORS[i])
            ax.plot(alpha_f_plus, np.abs(firstBornPluss)**2, ls="--", lw=0.3, color=cs.COLORS[i])
            ax.plot(alpha_f_minus, np.abs(firstBornMinus)**2, ls="--", lw=0.3, color=cs.COLORS[i])

        fname = "Figures/dwbaTotalPattern.pdf" 
        ax.set_xlabel("$\\alpha_f/\\alpha_c$")
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_yscale("log")
        ax.text(0.7, 0.17, "DWBA (solid)", transform=ax.transAxes)
        ax.text(0.7, 0.1, "BA (dashed)", transform=ax.transAxes)
        ax.legend(loc="lower left", frameon=False)
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname)) 
        

# Add a main function for plotting the total spectrum
def main():
    angles = np.array([0.5,1.0,2.0,10.0])
    try:
        combHandler = CombinedDetector()
        combHandler.plotTotal(angles)
    except Exception as e:
        print str(e)

if __name__ == "__main__":
    main()
