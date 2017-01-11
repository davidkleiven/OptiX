import tkinter as tk
import drawer as draw

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
        for i in range(self.showCount, len(self.wgs) ):
            self.drawer.addArc( float( self.wgs[i].rIn.get() ), float( self.wgs[i].angleIn.get() ),
                              self.wgs[i].curvature.get() )
            self.showCount += 1
        self.drawer.show()

    def save( self ):
        print ("Not implemented yet!")

    def pack( self ):
        for row in range(0, len(self.wgs)):
            self.wgs[row].grid(row)

        curRow = len(self.wgs)
        self.addButton.grid(row=curRow, column=0)
        self.removeButton.grid(row=curRow, column=1)
        self.showButton.grid(row=curRow, column=2)
        self.saveButton.grid(row=curRow, column=3)
