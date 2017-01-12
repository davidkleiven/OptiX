import tkinter as tk
import drawer as draw
import json
import numpy as np

class WgEntry:
    def __init__(self, master):
        self.master = master
        self.rLabel = tk.Label(self.master, text="R:")
        self.rIn = tk.Entry( self.master, width=10 )
        self.angleLabel = tk.Label( self.master, text="angle:")
        self.angleIn = tk.Entry( self.master, width=10 )
        self.curvature = tk.StringVar( self.master )
        self.availableCurvature = ["concave", "convex"]
        self.curvature.set( self.availableCurvature[0] )
        self.curvatureMenu = tk.OptionMenu( self.master, self.curvature, *self.availableCurvature )

        # Insert default values
        self.rIn.insert( tk.END, 40.0 )
        self.angleIn.insert( tk.END, 1.0 )

    def grid( self, row ):
        self.rLabel.grid(row=row, column=0)
        self.rIn.grid(row=row, column=1)
        self.angleLabel.grid(row=row, column=2)
        self.angleIn.grid(row=row, column=3)
        self.curvatureMenu.grid(row=row, column=4)

    def destroy( self ):
        self.rLabel.destroy()
        self.rIn.destroy()
        self.angleLabel.destroy()
        self.angleIn.destroy()

class Editor:
    def __init__( self, master ):
        self.frame = tk.Frame(master)
        self.frame.grid(row=0,column=0)
        self.drawer = draw.Drawer()

        self.wgs = []
        self.wgs.append( WgEntry(self.frame) )

        self.fnameEntry = tk.Entry( self.frame, width=10 )
        self.fnameEntryLabel = tk.Label( self.frame, text="Geofile" )
        self.fnameEntry.insert( tk.END, "InputFiles/waveguideGeometry%d.json"%( np.random.rand()*10000 ) )

        self.addButton = tk.Button( self.frame, text="AddWG", command=self.addWG )
        self.removeButton = tk.Button( self.frame, text="Remove last WG", command=self.removeLastWG )
        self.showButton = tk.Button( self.frame, text="Show", command=self.show )
        self.saveButton = tk.Button( self.frame, text="Save", command=self.save)
        self.pack()
        self.showCount = 0

    def addWG( self ):
        self.wgs.append( WgEntry(self.frame) )
        self.pack()

    def removeLastWG( self ):
        if ( len(self.wgs) == 1 ):
            return
        self.wgs[-1].destroy()
        self.wgs = self.wgs[:-1]
        self.pack()

    def show( self ):
        self.drawer.reset()
        for i in range(0, len(self.wgs) ):
            self.drawer.addArc( float( self.wgs[i].rIn.get() ), float( self.wgs[i].angleIn.get() ),
                              self.wgs[i].curvature.get() )
            self.showCount += 1
        self.drawer.show()

    def save( self ):
        fname = str( self.fnameEntry.get() )
        if ( fname[-4:] != "json" ):
            print ("File extension has to be .json!")
            return

        self.drawer.save( fname )
        geo = {}
        geo["figname"] = self.drawer.createFigname( fname )
        geo["geofile"] = fname
        geo["waveguides"] = []

        for wg in self.wgs:
            wgDict = {}
            wgDict["radius"] = float( wg.rIn.get() )
            wgDict["angle"] = float( wg.angleIn.get() )
            wgDict["curvature"] = wg.curvature.get()
            geo["waveguides"].append(wgDict)

        json_data = json.dumps( geo, indent=2, sort_keys=True )
        outfile = open( fname, 'w' )
        outfile.write( json_data )
        outfile.close()
        print ( "Geometry written to %s"%(fname) )

    def pack( self ):
        for row in range(0, len(self.wgs)):
            self.wgs[row].grid(row)

        curRow = len(self.wgs)
        self.addButton.grid(row=curRow, column=0)
        self.removeButton.grid(row=curRow, column=1)
        self.showButton.grid(row=curRow, column=2)
        self.saveButton.grid(row=curRow, column=3)
