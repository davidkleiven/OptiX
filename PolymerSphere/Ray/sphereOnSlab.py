import numpy as np
from matplotlib import pyplot as plt
import os
import sys
from xrt.backends.raycing import materials as rm
from xrt.backends.raycing import oes as roe
from xrt.backends import raycing
from xrt.backends.raycing import sources as rs
from xrt import plotter as xrtp
from xrt import runner as xrtr
from xrt.backends.raycing import run as rr
from xrt.backends.raycing import screens as rsc
from xrt.backends.raycing import oes

IS_PARAMETRIC = True
RADIUS = 1.0

# Using parametrised coordinates
class Sphere(roe.OE):
    def __init__(self, *args, **kwargs):
        self.Rm = RADIUS
        roe.OE.__init__(self,*args,**kwargs)

    def local_r(self, s, phi):
        return self.Rm

    def local_n(self, s, phi): 
        x,y,z = self.param_to_xyz(s,phi, self.Rm)
        n = np.zeros(3)
        nhat = np.array([x,y,z])/self.Rm
        return nhat

    def param_to_xyz(self, s, phi, r):
        x = r*np.sin(phi)*np.sin(s)
        y = r*np.cos(phi)*np.sin(s)
        z = r*np.cos(s)
        return x,y,z

    def xyz_to_param(self, x, y, z ):
        r = np.sqrt(x**2+y**2+z**2)
        s = np.arccos(z/r)
        phi = np.arccos(y/np.sqrt(x**2+y**2)) 
        return s,phi,r

def build_beamline(radius, nrays=raycing.nrays):
    beamline = raycing.BeamLine(height=0.0)
    beamline.fsm1 = rsc.Screen( beamline, "DiamondFSM1", (0.0,5.0*RADIUS,0.0))
    sigma = 5.0*radius
    E0 = 7300 # Energy in eV
    source = rs.GeometricSource(beamline, "GeometricSource", center=(0,-5.0*RADIUS,0),
             nrays=nrays, distx="normal", dx=sigma, distz="normal", dz=sigma, 
             distxprime="flat", dxprime=0.0, distzprime="flat", dzprime=0.0, 
             distE="lines", energies=(E0,), polarization="horizontal")
    fname = "sphere"
    limPhysX = [-radius, radius]
    limPhysY = [-radius, radius]
    material = rm.Material("Au", rho=19.3)
    beamline.sphere = Sphere(beamline, fname, [0.0,0.0,0.0], pitch=0.0, material=material, limPhysX=limPhysX, limPhysY=limPhysY,
    isParametric=IS_PARAMETRIC)

    beamline.fsm2 = rsc.Screen( beamline, "DiamondFSM2", (0.0, 1000.0*RADIUS, 0.0) )
    return beamline, fname

def run_process( beamline, runonlyfirst=False ):
    beamsource = beamline.sources[0].shine()
    beamFSM1 = beamline.fsm1.expose(beamsource)
    beamSphereGlobal, beamSphereLocalN = beamline.sphere.multiple_reflect( beamsource, maxReflections=100)
    beamFSM2 = beamline.fsm2.expose( beamSphereGlobal )
    outDict = {"beamSource":beamsource, "beamFSM1":beamFSM1, "beamSphereGlobal":beamSphereGlobal,
               "beamSphereLocalN":beamSphereLocalN, "beamFSM2":beamFSM2}
    return outDict

rr.run_process = run_process

def definePlots(beamline, fname):
    plots = []

    # First plot: Near field
    plot = xrtp.XYCPlot( "beamFSM1", (1,), xaxis=xrtp.XYCAxis(r"$x$", "$\mu m$"), yaxis=xrtp.XYCAxis(r"$z$", "$\mu m$" ), 
                         title="FSM1")
    plot.caxis.fwhmFormStr = None
    plot.saveName = ["beamFSM1.png",]
    plots.append(plot)

    # Scattered field:
    plot = xrtp.XYCPlotWithNumerOfReflections("beamFSM2", (1,), xaxis=xrtp.XYCAxis(r"$x$", "mm"), 
           yaxis=xrtp.XYCAxis(r"$z$", "mm"), caxis=xrtp.XYCAxis("Number of Reflections", "", bins=32, ppb=8,
           data=raycing.get_reflection_number), title="FSM2_Es")
    plot.saveName=["beamFSM2.png",]
    plots.append(plot)
   
    return plots  
    
def main():
    beamline, fname = build_beamline(RADIUS)
    plots = definePlots( beamline, fname )
    xrtr.run_ray_tracing(plots, repeats=40, updateEvery=1, beamLine=beamline, processes="half")
     
if __name__ == "__main__":
    main()
