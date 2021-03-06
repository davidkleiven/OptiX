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
        ax.text(0.7,0.1, "$\\alpha = %.1f\\alpha_c$"%(self.grazingAngle), transform=ax.transAxes)
        ax.legend(loc="lower left", ncol=2, frameon=False)
        fname = "Figures/termsTransmission.pdf"
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname)) 

    def adjustX(self, incAngle, scatteringAngles):
        alpha_f = self.detectorTransform.alpha_f(incAngle, scatteringAngles)
        alpha_fDeg = alpha_f*self.alpha_c*180.0/np.pi
        alpha_f = alpha_f[alpha_fDeg > -90.0]
        scatteringAngles = scatteringAngles[alpha_fDeg > -90.0]

        self.x = self.detectorPosition*np.tan(alpha_f*self.alpha_c)
        scatteringAngles = scatteringAngles[alpha_f<0.0]
        self.x = self.x[alpha_f<0.0]
        alpha_f = alpha_f[alpha_f<0.0]
        return scatteringAngles

    def plotPattern(self, angles):
        self.prepareDWBA(angles[0]) 
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for i in range(0, len(angles)): 
            if ( angles[i] < 1.0 ):
                label = "$\\alpha = %.1f\\alpha_c$"%(angles[i])
            elif ( angles[i] == 1 ):
                label = "$\\alpha = \\alpha_c$"%(angles[i])
            else:
                label = "$\\alpha=%d\\alpha_c$"%(angles[i]) 
            scatteringAngles = np.linspace(-89.9*np.pi/(180.0*self.alpha_c), 85.0*np.pi/(180.0*self.alpha_c), 100001)
            validScatAngles = self.adjustX(angles[i], scatteringAngles)
            self.prepareDWBA(angles[i]) 
            ba = self.born()
            tot = self.total()
            ax.plot(validScatAngles, np.abs(tot)**2, color=cs.COLORS[i], label=label)
            if ( i == (len(angles)-1) ):
                ax.plot(validScatAngles, np.abs(ba)**2, lw=0.3, color=cs.COLORS[len(angles)], label="BA")
            else:
                ax.plot(validScatAngles, np.abs(ba)**2, lw=0.3, color=cs.COLORS[len(angles)])

        fname = "Figures/dwbaTransmissionPattern.pdf"
        ax.set_xlabel(self.detectorTransform.axisLabel())
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1E-8)
        ax.legend(loc="upper left", frameon=False, labelspacing=0.05 )
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname)) 
