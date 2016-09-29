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
import transformDetector as td # Wrapper wich shifts the observation angle

class GrazingHandler:
    def __init__(self, usefilm):
        self.filmThickIsSet = False
        if ( usefilm ):
            self.coeff = scat.FilmScatterer()
        else:
            self.coeff = scat.PlaneScatterer()
            self.filmThickIsSet = True # Set to true since it is irrelevant
        self.detectorPosition = 100.0
        self.Rsphere = 1.0
        self.kR = 5.0
        self.usefilm = usefilm

        self.epsSubstIsSet = False
        self.prepDWBAIsCalled = False
        self.prepMSG = "You have to call prepareDWBA before doing any computations!"

        # Variables for storing computed quantities
        self.f1 = None # Sphere --> detector
        self.f2 = None # Substrate --> Sphere --> detector
        self.f3 = None # Sphere --> Substrate --> detector
        self.f4 = None # Substrate --> Sphere --> Substrate --> detector
        self.tot = None # Sum of all the above

        # Scattered wave vector
        self.kx_scat = None
        self.ky_scat = None
        self.kz_scat = None
        self.waveVector = np.zeros(3)
        self.q_parallel = None
        self.k = None # Think this is the same as kR so it is unnessecary, carries over from earlier development
        self.alpha_c = None
        self.x = np.linspace(0.0,2.5*self.detectorPosition, 100001)
        self.grazingAngle = 0.0

        self.detectorTransform = td.DetectorTransformer() # Default is the identity

    def setEpsilonSubst(self, epsSubst):
        self.coeff.eps2 = epsSubst
        self.epsSubstIsSet = True

    def setFilmThickness(self, filmThick):
        self.coeff.thickness = filmThick
        self.filmThickIsSet = True

    def prepareDWBA(self, grazingAngle): 
        if ( not self.filmThickIsSet ):
            raise Exception("You have to set the film thickness before computing DWBA quantities!")
        if ( not self.epsSubstIsSet ):
            raise Exception("You have to specify epsilon in the scatterer(s) before computing DWBA quantities!")
        mu = 1.0
        self.coeff.k = self.kR
        self.alpha_c = fe.criticalGrazingAngle( 1.0, np.sqrt(self.coeff.eps2*mu) )
        self.waveVector[0] = -self.kR*np.sin(grazingAngle*self.alpha_c)
        self.waveVector[2] = self.kR*np.cos(grazingAngle*self.alpha_c)
        self.grazingAngle = grazingAngle
        y = 0.0
        self.k = np.sqrt( np.sum(self.waveVector**2) )
        #self.q_parallel = np.sqrt( self.waveVector[1]**2 + self.waveVector[2]**2 ) 
        rHat_x = self.x/np.sqrt(self.x**2+y**2+self.detectorPosition**2)
        rHat_z = self.detectorPosition/np.sqrt(self.x**2+y**2+self.detectorPosition**2)
        self.kx_scat = self.k*rHat_x
        self.ky_scat = 0.0
        self.kz_scat = self.k*rHat_z
        self.q_parallel = np.sqrt( (self.ky_scat-self.waveVector[1])**2 + (self.kz_scat-self.waveVector[2])**2 ) 
        self.prepDWBAIsCalled = True

    def born(self):
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)

        qx = self.kx_scat - self.waveVector[0]
        qR = np.sqrt( self.q_parallel**2 + qx**2 )
        qR[qx<0.0] = -qR[qx<0.0]
        self.f1 = self.formFactorSphere(qR)
        return self.f1

    def bornSubstrateSphere(self): 
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)
        qx = self.kx_scat + self.waveVector[0]
        qR = np.sqrt( self.q_parallel**2 + qx**2 )
        angleWithSubstrateDeg = np.arccos(-self.waveVector[0]/self.k)*180.0/np.pi
        self.f2 = self.formFactorSphere(qR)*self.coeff.reflection(angleWithSubstrateDeg)
        return self.f2
    
    def bornSphereSubstrate(self): 
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)
        qx = -self.kx_scat - self.waveVector[0]
        qR = np.sqrt( self.q_parallel**2 + qx**2 )
        angleWithSubstrateDeg = np.arccos(self.kx_scat/self.k)*180.0/np.pi
        self.f3 = self.formFactorSphere(qR)*self.coeff.reflection(angleWithSubstrateDeg)
        self.f3[-self.kx_scat>0.0] = 0.0 # Ray leaves sphere with kx = -kx_scat
        return self.f3

    def bornSubstrateSphereSubstrate(self): 
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)
        qx = -self.kx_scat + self.waveVector[0]
        qR = np.sqrt( self.q_parallel**2 + qx**2 )
        angleWithSubstrateFirstDeg = np.arccos(np.abs(self.waveVector[0])/self.k)*180.0/np.pi
        angleWithSubstrateSecondDeg = np.arccos(self.kx_scat/self.k)*180.0/np.pi
        self.f4 = self.formFactorSphere(qR)*self.coeff.reflection(angleWithSubstrateFirstDeg)*self.coeff.reflection(angleWithSubstrateSecondDeg)
        self.f4[-self.kx_scat>0.0] = 0.0
        return self.f4

    def bornTotal(self):
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)
        self.born()
        self.bornSubstrateSphere()
        self.bornSphereSubstrate()
        self.bornSubstrateSphereSubstrate()
        return self.f1 + self.f2 + self.f3 + self.f4

    def formFactorSphere(self, qR ):
        form = ( np.sin(qR) - qR*np.cos(qR) )/qR**3
        return form

    def plotTerms(self): 
        if ( not self.prepDWBAIsCalled ):
            raise Exception(self.prepMSG)
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
            ax.text( 0.1, 0.15, "$\\alpha_i = %d\\alpha_c$"%(self.grazingAngle), transform=ax.transAxes)
        elif( self.grazingAngle == 1 ):
            ax.text( 0.1,0.15, "$\\alpha_i = \\alpha_c$", transform=ax.transAxes)
        else:
            ax.text( 0.1, 0.15, "$\\alpha_i = %.1f\\alpha_c$"%(self.grazingAngle), transform=ax.transAxes)
        ax.legend( loc="lower right", frameon=False, ncol=4 )
        
        if ( self.usefilm ):
            fname = "Figures/patternFilm1D.pdf"
        else:
            fname = "Figures/pattern1D.pdf"

        fig.savefig( fname, bbox_inches="tight" )
        print ("Figure written to %s"%(fname))

    def totalAngleSweep(self, angles):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for i in range(0, len(angles)):
            if ( angles[i] < 1.0 ):
                label = "$\\alpha = %.1f\\alpha_c$"%(angles[i])
            elif ( angles[i] == 1 ):
                label = "$\\alpha = \\alpha_c$"%(angles[i])
            else:
                label = "$\\alpha=%d\\alpha_c$"%(angles[i])
            self.prepareDWBA(angles[i])
            alpha_f = np.arctan(self.x/self.detectorPosition)
            alpha_f /= self.alpha_c
            scatteringAngle = self.detectorTransform.scatteringAngle(angles[i], alpha_f)
            tot = np.abs(self.bornTotal()**2)
            firstBorn = np.abs(self.f1)**2
            ax.plot( scatteringAngle, tot, color=cs.COLORS[i], label=label)
            if ( i == (len(angles)-1) ):
                ax.plot( scatteringAngle, firstBorn, lw=0.3, color=cs.COLORS[len(angles)], label="BA")
            else:
                ax.plot( scatteringAngle, firstBorn, lw=0.3, color=cs.COLORS[len(angles)])
        ax.set_xlabel( self.detectorTransform.axisLabel() )
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1E-8)
        ax.legend(loc="upper right", frameon=False, labelspacing=0.2)
        if ( self.usefilm ):
            fname = "Figures/dwbaPatternFilm.pdf"
        else:
            fname = "Figures/dwbaPattern.pdf"
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))
