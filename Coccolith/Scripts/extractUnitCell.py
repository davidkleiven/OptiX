import sys
import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
#from matplotlib import pyplot as plt
from matplotlib import figure as mfig
from matplotlib.backends import backend_tkagg as btk
import tkinter as tk
from matplotlib import patches
from tkinter import filedialog
from scipy.ndimage import interpolation as imginter

class Application(tk.Tk):
    def __init__( self, *args, **kwargs ):
        tk.Tk.__init__(self,*args,**kwargs)
        tk.Tk.wm_title(self,"Unit Cell Extractor")

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)

        self.frames = {}
        self.frames[PlotFrame] = PlotFrame(container,self)
        self.frames[PlotFrame].grid(row=0,column=0)
        self.frames[PlotFrame].tkraise()

class PlotFrame( tk.Frame ):
    def __init__( self, parent, controller ):
        tk.Frame.__init__( self, parent )
        self.data = []

        self.fig1 = mfig.Figure(figsize=(8,8), dpi=100)
        self.ax1 = self.fig1.add_subplot(1,2,1)
        self.ax2 = self.fig1.add_subplot(1,2,2)
        self.canvas = btk.FigureCanvasTkAgg(self.fig1, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = btk.NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Add entry fields
        subframe = tk.Frame(self)
        subframe.pack(side=tk.TOP)
        tk.Label(subframe, text="xmin").grid(row=0,column=0)
        self.xmin = tk.Entry(subframe)
        self.xmin.grid(row=0,column=1)
        self.xmin.insert(0,0)
        tk.Label(subframe, text="xmax").grid(row=0,column=2)
        self.xmax = tk.Entry(subframe)
        self.xmax.grid(row=0,column=3)
        self.xmax.insert(0,1)

        tk.Label(subframe, text="ymin").grid(row=1,column=0)
        self.ymin = tk.Entry(subframe)
        self.ymin.grid(row=1,column=1)
        self.ymin.insert(0,0)

        tk.Label(subframe, text="ymax").grid(row=1,column=2)
        self.ymax = tk.Entry(subframe)
        self.ymax.grid(row=1, column=3)
        self.ymax.insert(0,1)

        tk.Label(subframe, text="rotation").grid(row=2,column=0)
        self.angle = tk.Entry(subframe)
        self.angle.grid(row=2,column=1)
        self.angle.insert(0,0)
        self.isFirst = True
        tk.Button(subframe, text="Expand Unit Cell", command=self.expandUnitCell).grid(row=3,column=1)
        tk.Button(subframe, text="Save Unit Cell", command=self.saveUnit).grid(row=3,column=2)
        tk.Button(subframe, text="rotate", command=self.rotateImag).grid(row=3,column=0)

        self.selectedXmin = 0
        self.selectedXmax = 1
        self.selectedYmin = 0
        self.selectedYmax = 1

        self.oldfname = "defaultFilename_100.raw"

    def expandUnitCell( self ):
        self.selectedXmin = int(self.xmin.get())
        self.selectedXmax = int(self.xmax.get())
        self.selectedYmin = int(self.ymin.get())
        self.selectedYmax = int(self.ymax.get())
        unitCell = self.getUnitCell()
        self.plotSuperCell( unitCell )
        self.plot()

    def plotSuperCell( self, unitCell, N=3 ):
        array = np.zeros( (unitCell.shape[0]*N, unitCell.shape[1]*N) )
        for i in range(0,N):
            for j in range(0,N):
                array[i*unitCell.shape[0]:(i+1)*unitCell.shape[0],j*unitCell.shape[1]:(j+1)*unitCell.shape[1]] = unitCell
        self.ax2.clear()
        self.ax2.imshow(array, cmap="bone", aspect="equal")
        self.canvas.draw()

    def getUnitCell( self ):
        ucell = self.data[self.selectedYmin:self.selectedYmax,self.selectedXmin:self.selectedXmax]
        return ucell

    def plot(self):
        if ( len(self.data) == 0 ):
            raise Exception("No data given!")
        self.ax1.clear()
        self.ax1.imshow(self.data, cmap="bone", aspect="equal")
        w = self.selectedXmax-self.selectedXmin
        h = self.selectedYmax-self.selectedYmin
        self.ax1.add_patch( patches.Rectangle((self.selectedXmin,self.selectedYmin), w, h, fill=False))
        if ( self.isFirst ):
            self.ax2.clear()
            self.ax2.imshow(self.data,cmap="bone", aspect="equal")
            self.isFirst = False
        self.canvas.draw()
        self.canvas.show()

    def saveUnit( self ):
        indx = self.oldfname.find("_")
        fname = self.oldfname[:indx]+"UnitCell"
        resolution = int(self.oldfname.split("_")[1])
        unitcell = self.getUnitCell()
        fname += "_%d_%d_%d.raw"%(resolution, unitcell.shape[0], unitcell.shape[1])
        unitcell.ravel().tofile(fname)
        print ("File saved in %s"%(fname))

    def rotateImag( self ):
        self.data = imginter.rotate( self.data, int(self.angle.get()), reshape=False )
        self.plot()
        self.expandUnitCell()
        self.angle.insert(0,0)


def getSize( fname ):
    splitted = fname.split(".")[0].split("_")
    if ( len(splitted) != 5 ):
        raise (Exception("Filename in wrong format: Expected <somename_resolution_Nx_Ny_Nz.raw") )

    res = int(splitted[-4])
    Nz = int(splitted[-1])
    Ny = int(splitted[-2])
    Nx = int(splitted[-3])
    return res,Nx,Ny,Nz

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python extractUnitCell.py <fname_res_Nx_Ny_Nz.raw>")
        return 1
    data = np.fromfile(argv[0], dtype=np.uint8)
    try:
        res,Nx,Ny,Nz = getSize(argv[0])
        data = data.reshape((Nx,Ny,Nz), order="F")
        app = Application()
        proj = data.sum(axis=1).T
        proj = proj*255/int(proj.max())
        app.frames[PlotFrame].data = proj.astype(np.uint8)
        app.frames[PlotFrame].oldfname = argv[0]
        del data
        app.frames[PlotFrame].plot()
        app.mainloop()
    except Exception as exc:
        print (str(exc))
        return 1
    return 0

if __name__ == "__main__":
    main( sys.argv[1:])
