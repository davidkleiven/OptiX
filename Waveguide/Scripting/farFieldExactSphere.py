import sys
import proj3D
import h5py as h5
from PLOD import controlGUI as cg
import numpy as np
import tkinter as tk
from scipy import special as spec
from scipy import integrate
import json

k = 40.01 # nm^{-1} wavenumber
#fname = "data/sphere828475.h5"

#fname = "data/sphere558434.h5"
fname = "data/sphere334411.h5"
fname = "data/sphere912348.h5"
fname = "data/sphere799650.h5"
fname = "data/sphere483220.h5"
fname = "data/sphere450350.h5"
fname = "data/sphere683666.h5"
fname = "data/sphere558434.h5"
fname = "data/sphere718369.h5"

class LayeredSphere:
    def __init__( self ):
        self.radii = []
        self.delta = []

    def formsphere( self, q, R ):
        V = (4.0*np.pi*R**3)/3.0
        f =  3.0*V*(np.sin(q*R) - q*R*np.cos(q*R))/(q*R)**3
        f[np.abs(q*R) < 1E-3] = V
        return f

    def formFactor( self, q ):
        assert( len(self.radii) == len(self.delta) )
        tot = self.delta[0]*self.formsphere( q, self.radii[0] )*self.delta[0]
        for i in range( 1, len(self.radii) ):
            tot += self.delta[i]*(self.formsphere( q, self.radii[i]) - self.formsphere( q, self.radii[i-1]) )
        return tot

def integrandReal( r, q, R, delta ):
    return r*spec.jn( 0, q*r )*( np.cos(-delta*k*np.sqrt(R**2 - r**2) ) - 1.0 )

def integrandImag( r, q, R, delta ):
    return r*spec.jn( 0, q*r )*np.sin(-delta*k*np.sqrt(R**2 - r**2) )

def formFactorCircularStop( q, R ):
    q[np.abs(q) < 1E-15] = 1E-15
    return 2.0*spec.j1( q*R )/(q*R )

def main( argv ):
    paramFile = open( argv[0], 'r' )
    params = json.load( paramFile )
    paramFile.close()

    plotCircularStop = 0
    try:
        plotCircularStop = params["plotCircularStop"]
    except:
        plotCircularStop = 0

    lim = proj3D.Limits()
    ffPlot = proj3D.FarField()
    with h5.File(params["datafile"], 'r') as hf:
        ff = hf.get("/data/farField")
        lim.qmin = ff.attrs.get("qmin")
        lim.qmax = ff.attrs.get("qmax")
        ffPlot.data = np.array( ff )

    ffPlot.limits = lim
    ffPlot.centerOfMass()

    root = tk.Tk()
    control = cg.Control( root )
    control = ffPlot.projectionX( control, color="#0571b0", label="Num" )
    ax = control.plots.axes[0].ax
    amp = np.max( ffPlot.data )

    q = np.linspace( lim.qmin, lim.qmax, 1001 )

    spheres = LayeredSphere()
    for entry in params["spheres"]:
        spheres.radii.append( entry["radius"] )
        spheres.delta.append( entry["delta"] )

    F = spheres.formFactor( q )**2

    F *= ( amp/np.max(F) )

    R = params["spheres"][0]["radius"]
    delta = params["spheres"][0]["radius"]

    Fc = formFactorCircularStop( q, R )**2
    Fc *= ( amp/np.max(Fc) )

    Fphase = np.zeros(len(q)) + 1j*np.zeros(len(q))
    for i in range(0,len(q)):
        Fphase[i] = integrate.quad( integrandReal, 0.0, R, args=(q[i],R,delta))[0] + 1j*integrate.quad( integrandImag, 0.0, R, args=(q[i],R,delta))[0]


    control.plots.axes[0].ax.plot( q, F, color="#ca0020", label=r"\$F(q)\$" )
    control.plots.axes[0].ax.legend( loc="upper left", frameon=False, labelspacing=0.01, borderpad=0.0, handletextpad=0.3, handlelength=0.1 )

    if ( plotCircularStop == 1 ):
        control.plots.axes[0].ax.plot( q, Fc, color="#4daf4a", label=r"\$F_c(q)\$")

    control.plots.axes[0].ax.legend( loc="upper left", frameon=False, labelspacing=0.01, borderpad=0.0, handletextpad=0.3, handlelength=0.1 )
    root.mainloop()

if __name__ == "__main__":
    main( sys.argv[1:] )
