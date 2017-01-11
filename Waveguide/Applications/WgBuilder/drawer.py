import tkinter as tk
import numpy as np

class Arc:
    def __init__(self):
        self.x0 = 0.0
        self.y0 = 0.0
        self.R = 1.0
        self.angle = 1.0

    def endPoint( self ):
        x1 = self.R*np.sin( self.angle*np.pi/180.0 )
        y1 = self.R-self.R*np.cos( self.angle*np.pi/180.0 )
        return x1, y1

class Drawer:
    def __init__( self, master ):
        self.frame = tk.Frame(master)
        self.frame.grid(row=0,column=1)
        self.dim=400
        self.extentX = 0.0
        self.extentY = 0.0
        self.initRadius = 1.0
        self.canvas = tk.Canvas( self.frame, width=self.dim, height=self.dim, bg="black" )
        self.canvas.pack()
        self.arcs = []
        self.sqSize = 2.0
        self.colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"]

        # Arc layout
        self.style = "arc"
        self.outline = "white"

    def crdToPixX( self, x ):
        return x*self.dim/self.extentX

    def crdToPixY( self, y ):
        return y*self.dim/self.extentY

    def addArc( self, R, angle ):
        arc = Arc()
        arc.R = R
        arc.angle = angle
        if ( len(self.arcs) == 0 ):
            arc.x0 = 0.0
            arc.y0 = 0.0
        else:
            arc.x0, arc.y0 = self.arcs[-1].endPoint()
            arc.x0 += self.arcs[-1].x0
            arc.y0 += self.arcs[-1].y0

        self.arcs.append(arc)

    def setExtent( self ):
        self.extentX = 0
        self.extentY = 0
        for arc in self.arcs:
            x1, y1 = arc.endPoint()
            if ( arc.x0+x1 > self.extentX ):
                self.extentX = arc.x0+x1

            if ( arc.y0+y1 > self.extentY ):
                self.extentY = arc.y0+y1

    def draw( self ):
        self.setExtent()
        self.canvas.create_rectangle( 0, 0, self.dim, self.dim, fill="black")
        counter = 0
        for arc in self.arcs:
            x1, y1 = arc.endPoint()
            x0 = self.crdToPixX( arc.x0 )
            y0 = self.crdToPixY( arc.y0 )
            x1 = self.crdToPixX( arc.x0+x1 )
            y1 = self.crdToPixY( arc.y0+y1 )

            print (x0, y0, x1, y1)
            self.canvas.create_line(x0, y0, x1, y1, fill=self.colors[counter%len(self.colors)], width=6)
            counter += 1
