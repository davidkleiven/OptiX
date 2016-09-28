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

class GrazingHandler:
    def __init__(self, usefilm):
        if ( usefilm ):
            self.coeff = scat.FilmScatterer()
        else:
            self.coeff = scat.PlaneScatterer()
        self.detectorPosition = 100.0
        self.Rsphere = 1.0
        self.kR = 5.0
        self.usefilm = usefilm

        # Variables for storing computed quantities
        self.f1 = None # Sphere --> detector
        self.f2 = None # Substrate --> Sphere --> detector
        self.f3 = None # Sphere --> Substrate --> detector
        self.f4 = None # Substrate --> Sphere --> Substrate --> detector
        self.tot = None # Sum of all the above

        # Scattered wave vector
        self.kx_scat = None
        self.waveVector = np.zeros(3)
        self.q_parallel = None
        self.k = None # Think this is the same as kR so it is unnessecary, carries over from earlier development
        self.alpha_c = None
        self.x = np.linspace(0.0,2.5*self.detectorPosition, 100001)
        self.grazingAngle = 0.0

    def setEpsilonSubst(self, epsSubst):
        self.coeff.eps2 = epsSubst
    def setFilmThickness(self, filmThick):
        self.coeff.thickness = filmThick

    def prepareDWBA(self, grazingAngle): 
        mu = 1.0
        self.coeff.k = self.kR
        self.alpha_c = fe.criticalGrazingAngle( 1.0, np.sqrt(self.coeff.eps2*mu) )
        self.waveVector[0] = -self.kR*np.sin(grazingAngle*self.alpha_c)
        self.waveVector[2] = self.kR*np.cos(grazingAngle*self.alpha_c)
        self.grazingAngle = grazingAngle
        y = 0.0
        self.k = np.sqrt( np.sum(self.waveVector**2) )
        self.q_parallel = np.sqrt( self.waveVector[1]**2 + self.waveVector[2]**2 ) 
        rHat_x = self.x/np.sqrt(self.x**2+y**2+self.detectorPosition**2)
        self.kx_scat = self.k*rHat_x

    def born(self):
        qx = self.kx_scat - self.waveVector[0]
        qR = np.sqrt( self.q_parallel**2 + qx**2 )
        qR[qx<0.0] = -qR[qx<0.0]
        self.f1 = self.formFactorSphere(qR)
        return self.f1

    def bornSubstrateSphere(self): 
        qx = self.kx_scat + self.waveVector[0]
        qR = np.sqrt( self.q_parallel**2 + qx**2 )
        angleWithSubstrateDeg = np.arccos(-self.waveVector[0]/self.k)*180.0/np.pi
        self.f2 = self.formFactorSphere(qR)*self.coeff.reflection(angleWithSubstrateDeg)
        return self.f2
    
    def bornSphereSubstrate(self): 
        qx = -self.kx_scat - self.waveVector[0]
        qR = np.sqrt( self.q_parallel**2 + qx**2 )
        angleWithSubstrateDeg = np.arccos(self.kx_scat/self.k)*180.0/np.pi
        self.f3 = self.formFactorSphere(qR)*self.coeff.reflection(angleWithSubstrateDeg)
        self.f3[-self.kx_scat>0.0] = 0.0 # Ray leaves sphere with kx = -kx_scat
        return self.f3

    def bornSubstrateSphereSubstrate(self): 
        qx = -self.kx_scat + self.waveVector[0]
        qR = np.sqrt( self.q_parallel**2 + qx**2 )
        angleWithSubstrateFirstDeg = np.arccos(np.abs(self.waveVector[0])/self.k)*180.0/np.pi
        angleWithSubstrateSecondDeg = np.arccos(self.kx_scat/self.k)*180.0/np.pi
        self.f4 = self.formFactorSphere(qR)*self.coeff.reflection(angleWithSubstrateFirstDeg)*self.coeff.reflection(angleWithSubstrateSecondDeg)
        self.f4[-self.kx_scat>0.0] = 0.0
        return self.f4

    def bornTotal(self):
        self.born()
        self.bornSubstrateSphere()
        self.bornSphereSubstrate()
        self.bornSubstrateSphereSubstrate()
        return self.f1 + self.f2 + self.f3 + self.f4

    def formFactorSphere(self, qR ):
        form = ( np.sin(qR) - qR*np.cos(qR) )/qR**3
        return form

    def plotTerms(self): 
        self.born()
        self.bornSubstrateSphere()
        self.bornSphereSubstrate()
        self.bornSubstrateSphereSubstrate()
        alpha_f = np.arctan(self.x/self.detectorPosition)/self.alpha_c
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( alpha_f, np.abs(self.f1)**2, color=cs.COLORS[0], label="A")
        ax.plot( alpha_f, np.abs(self.f2)**2, color=cs.COLORS[1], label="B")
        ax.plot( alpha_f, np.abs(self.f3)**2, color=cs.COLORS[2], label="C")
        ax.plot( alpha_f, np.abs(self.f4)**2, color=cs.COLORS[3], label="D")
        ax.set_yscale('log')
        ax.set_xlabel("$\\alpha_f/\\alpha_c$")
        ax.set_ylabel("Intensity (a.u.)")
        if ( self.grazingAngle > 1.0 ):
            ax.text( 0.8*alpha_f[-1], 1E-5, "$\\alpha_i = %d\\alpha_c$"%(self.grazingAngle))
        elif( self.grazingAngle == 1 ):
            ax.text( 0.8*alpha_f[-1], 1E-5, "$\\alpha_i = \\alpha_c$") 
        else:
            ax.text( 0.8*alpha_f[-1], 1E-5, "$\\alpha_i = %.1f\\alpha_c$"%(self.grazingAngle))
        ax.legend( loc="lower right", frameon=False, ncol=4 )
        
        if ( self.usefilm ):
            fname = "Figures/patternFilm1D.pdf"
        else:
            fname = "Figures/pattern1D.pdf"

        fig.savefig( fname, bbox_inches="tight" )
        print ("Figure written to %s"%(fname))
