import sys
sys.path.append("../../FresnelBEM")
import numpy as np
import fresnelExact as fe
import filmCoefficients as fc

class PlaneScatterer:
    '''
    @brief Class for wrapper scattering from plane
    '''
    def __init__(self):
        self.eps1 = 1.0
        self.eps2 = 1.0
        self.mu1 = 1.0
        self.mu2 = 1.0
        self.k = 1.0
        self.polarisation = "TE"

    def reflection(self, incAngleDeg):
        if ( self.polarisation == "TE" ):
            return fe.rs(self.eps1, self.eps2, self.mu1, self.mu2, incAngleDeg, self.k) 
        return fe.rp(self.eps1, self.eps2, self.mu1, self.mu2, incAngleDeg, self.k) 

    def transmission(self, incAngleDeg):
        if ( self.polarisation == "TE" ):
            return fe.ts(self.eps1, self.eps2, self.mu1, self.mu2, incAngleDeg, self.k) 
        return fe.tp(self.eps1, self.eps2, self.mu1, self.mu2, incAngleDeg, self.k) 

class FilmScatterer(PlaneScatterer):
    '''
    @brief Wrapper for the film coefficients
    '''
    def __init__(self):
        PlaneScatterer.__init__(self)
        self.thickness = 1.0
    
    def bothCoeff(self, incAngleDeg):
        n1 = np.sqrt(self.eps1*self.mu1)
        n2 = np.sqrt(self.eps2*self.mu2)
        incGrazingAngle = 90.0-incAngleDeg
        critGrazingAngle = fe.criticalGrazingAngle(n1,n2)
        incAngleInUnitsOfAlpha_c = incGrazingAngle*np.pi/(180.0*critGrazingAngle)
        return fc.filmCoefficients(n1,n2, self.thickness, self.k, incAngleInUnitsOfAlpha_c, pol=self.polarisation)

    def reflection(self, incAngle):
        r,t = self.bothCoeff(incAngle)
        return r
    def transmission(self, incAngleDeg):
        r,t = self.bothCoeff(incAngle)
        return t
            
        
