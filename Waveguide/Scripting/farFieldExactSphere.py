import sys
import proj3D
import h5py as h5
from PLOD import controlGUI as cg
import numpy as np
import tkinter as tk
from scipy import special as spec
from scipy import integrate

R = 2000.0 # nm
k = 40.01 # nm^{-1} wavenumber
delta = 8.90652E-6
#fname = "data/sphere828475.h5"
fname = "data/sphere558434.h5"
def formsphere( q ):
    f =  (np.sin(q*R) - q*R*np.cos(q*R))/(q*R)**3
    f[np.abs(q*R) < 1E-3] = 1.0/3.0
    return f

def integrandReal( r, q ):
    return r*spec.jn( 0, q*r )*( np.cos(-delta*k*np.sqrt(R**2 - r**2) ) - 1.0 )

def integrandImag( r, q ):
    return r*spec.jn( 0, q*r )*np.sin(-delta*k*np.sqrt(R**2 - r**2) )

def formFactorCircularStop( q ):
    return 2.0*spec.jn( 1, q*R )/(q*R )

def main( argv ):
    plotCircularStop = 0
    if ( len(argv) > 0 ):
        plotCircularStop = argv[0]

    lim = proj3D.Limits()
    ffPlot = proj3D.FarField()
    with h5.File(fname, 'r') as hf:
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
    F = formsphere( q )**2
    F *= ( amp/np.max(F) )

    Fc = formFactorCircularStop( q )**2
    Fc *= ( amp/np.max(Fc) )

    Fphase = np.zeros(len(q)) + 1j*np.zeros(len(q))
    for i in range(0,len(q)):
        Fphase[i] = integrate.quad( integrandReal, 0.0, R, args=(q[i],))[0] + 1j*integrate.quad( integrandImag, 0.0, R, args=(q[i],))[0]


    control.plots.axes[0].ax.plot( q, F, color="#ca0020", label=r"\$F(q)\$")
    control.plots.axes[0].ax.legend( loc="upper left", frameon=False, labelspacing=0.01, borderpad=0.0, handletextpad=0.3, handlelength=0.1 )

    if ( plotCircularStop == 1 ):
        control.plots.axes[0].ax.plot( q, Fc, color="#4daf4a", label=r"\$F_c(q)\$")
    control.plots.axes[0].ax.legend( loc="upper left", frameon=False, labelspacing=0.01, borderpad=0.0, handletextpad=0.3, handlelength=0.1 )
    root.mainloop()

if __name__ == "__main__":
    main( sys.argv[1:] )
