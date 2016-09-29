import grazing as graz
import sys
sys.path.append("../../FresnelBEM")
sys.path.append("../../FresnelFDTD")
sys.path.append("../../")
import matplotlib
import mplLaTeX as ml
matplotlib.rcParams.update( ml.params )
import numpy as np
from matplotlib import pyplot as plt
import json
import fresnelExact as fe
import colorScheme as cs
import scatteringStructures as scat


class GrazingTransmissionHandler(graz.GrazingHandler):
    def __init__(self):
        graz.GrazingHandler.__init__(self, True) # usefilm True. Transmission does not make sense for infinite thick plane

        # Flip x as we are now looking below the film
        self.x = -self.x
        
        # Computed quantities
        self.ba = None
        self.f1 = None
        self.f2 = None
        
    def bornSphereSubstrate(self):
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)
        self.ba = self.born() 
        angleWithSubstrateDeg = np.arccos(self.kx_scat/self.k)*180.0/np.pi
        self.f1 = self.ba*self.coeff.transmission(angleWithSubstrateDeg)
        return self.f1

    def bornSubstrateSphere(self):
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)
        substSph = graz.GrazingHandler.bornSubstrateSphere(self)
        angleWithSubstrateDeg = np.arccos(self.kx_scat/self.k)*180.0/np.pi
        self.f2 = substSph*self.coeff.transmission(angleWithSubstrateDeg)
        return self.f2

    def total(self):
        self.bornSphereSubstrate()
        self.bornSubstrateSphere()
        return self.f1 + self.f2 

    def plotTerms(self):
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)
        self.bornSphereSubstrate()
        self.bornSubstrateSphere() 
        alpha_f = np.arctan(self.x/self.detectorPosition)/self.alpha_c

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(alpha_f, np.abs(self.f1)**2, color=cs.COLORS[0], label="A")  
        ax.plot(alpha_f, np.abs(self.f2)**2, color=cs.COLORS[4], label="B")
        ax.set_xlabel("$\\alpha_f/\\alpha_c$")
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_yscale("log")
        ax.legend(loc="lower left", ncol=2, frameon=False)
        fname = "Figures/termsTransmission.pdf"
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname)) 

    def plotPattern(self, angles):
        self.prepareDWBA(angles[0]) 
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        alpha_f = np.arctan(self.x/self.detectorPosition)/self.alpha_c
        for i in range(0, len(angles)): 
            if ( angles[i] < 1.0 ):
                label = "$\\alpha = %.1f\\alpha_c$"%(angles[i])
            elif ( angles[i] == 1 ):
                label = "$\\alpha = \\alpha_c$"%(angles[i])
            else:
                label = "$\\alpha=%d\\alpha_c$"%(angles[i])
            
            self.prepareDWBA(angles[i]) 
            ba = self.born()
            tot = self.total()
            ax.plot(alpha_f, np.abs(ba)**2, ls="--", lw=0.3, color=cs.COLORS[i])
            ax.plot(alpha_f, np.abs(tot)**2, color=cs.COLORS[i], label=label)
        fname = "Figures/dwbaTransmissionPattern.pdf"
        ax.set_xlabel("$\\alpha_f/\\alpha_c$")
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_yscale("log")
        ax.text(0.7, 0.17, "DWBA (solid)", transform=ax.transAxes)
        ax.text(0.7, 0.1, "BA (dashed)", transform=ax.transAxes)
        ax.legend(loc="lower left", frameon=False)
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname)) 
