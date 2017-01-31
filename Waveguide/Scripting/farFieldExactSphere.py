import proj3D
import h5py as h5
from PLOD import controlGUI as cg
import numpy as np
import tkinter as tk

R = 2000.0 # nm
fname = "data/sphere828475.h5"

def formsphere( q ):
    f =  (np.sin(q*R) - q*R*np.cos(q*R))/(q*R)**3
    f[np.abs(q*R) < 1E-3] = 1.0/3.0
    return f

def main():
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
    control = ffPlot.projectionX( control )
    ax = control.plots.axes[0].ax
    amp = np.max( ffPlot.data )

    q = np.linspace( lim.qmin, lim.qmax, 1001 )
    F = formsphere( q )**2
    F *= ( amp/np.max(F) )
    print ( np.max(F))
    control.plots.axes[0].ax.plot( q*10.0, F, color="#377eb8")
    root.mainloop()

if __name__ == "__main__":
    main()
