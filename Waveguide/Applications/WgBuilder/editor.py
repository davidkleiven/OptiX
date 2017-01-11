import tkinter as tk

class WgEntry:
    def __init__(self, master):
        self.master = master
        self.rLabel = tk.Label(self.master, text="R:")
        self.rIn = tk.Entry( self.master, width=10 )
        self.angleLabel = tk.Label( self.master, text="angle:")
        self.angleIn = tk.Entry( self.master, width=10 )

        # Insert default values
        self.rIn.insert( tk.END, 40.0 )
        self.angleIn.insert( tk.END, 1.0 )

    def grid( self, row ):
        self.rLabel.grid(row=row, column=0)
        self.rIn.grid(row=row, column=1)
        self.angleLabel.grid(row=row, column=2)
        self.angleIn.grid(row=row, column=3)

    def destroy( self ):
        self.rLabel.destroy()
        self.rIn.destroy()
        self.angleLabel.destroy()
        self.angleIn.destroy()

class Editor(tk.Frame):
    def __init__( self, master ):
        tk.Frame.__init__(self, master)

        self.wgs = []
        self.wgs.append( WgEntry(self.master) )

        self.addButton = tk.Button( self.master, text="AddWG", command=self.addWG )
        self.removeButton = tk.Button( self.master, text="Remove last WG", command=self.removeLastWG )
        self.showButton = tk.Button( self.master, text="Show", command=self.show )
        self.saveButton = tk.Button( self.master, text="Save", command=self.save)
        self.pack()

    def addWG( self ):
        self.wgs.append( WgEntry(self.master) )
        self.pack()

    def removeLastWG( self ):
        if ( len(self.wgs) == 1 ):
            return
        self.wgs[-1].destroy()
        self.wgs = self.wgs[:-1]
        self.pack()

    def show( self ):
        print ("Not implemented yet!")

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
        self.grid(row=0, column=0)
