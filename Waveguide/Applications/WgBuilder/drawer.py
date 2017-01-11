import tkinter as tk
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches as ptc

class Arc:
    def __init__(self):
        self.x0 = 0.0
        self.y0 = 0.0
        self.R = 1.0
        self.angle = 1.0

    def endPoint( self ):
        x1 = self.R*np.sin( self.angle*np.pi/180.0 )
        y1 = self.R*( 1.0-np.cos( self.angle*np.pi/180.0 ) )
        y1 = self.R*2.0*np.sin(0.5*self.angle*np.pi/180.0)**2
        return x1, y1

class Drawer:
    def __init__( self ):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1,1,1)
        self.x0 = 0.0
        self.y0 = 0.0
        self.xc = 0.0
        self.yc = 0.0
        self.xmin = 0.0
        self.xmax = 0.0
        self.ymin = 0.0
        self.ymax = 0.0
        self.tanAlpha = 0.0
        self.angle = 90.0
        self.ax.plot([0,0],[0,0])
        self.isFirst = True
        self.colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"]
        self.current = 0
        self.currentCurvature = "concave"
        self.ax.set_axis_bgcolor( "#404040" )
        self.fig.canvas.set_window_title("Wave Guide Geometry ")
        self.fig.patch.set_visible(False)
        self.annotations = []


    def customizePlot( self ):
        self.ax.spines["left"].set_position("zero")
        self.ax.spines["right"].set_color("none")
        self.ax.spines["left"].set_smart_bounds(True)
        self.ax.spines["bottom"].set_smart_bounds(True)
        self.ax.yaxis.set_ticks_position("left")
        self.ax.xaxis.set_ticks_position("top")
        self.ax.yaxis.set_ticks([])
        self.ax.xaxis.set_ticks([])

        dx = self.xmax - self.xmin
        dy = self.ymax - self.ymin
        ypos = self.ymax - 0.1*dy
        xpos = self.xmin + 0.05*dx
        self.annotations.append( self.ax.annotate(s="%.1f mm"%(dx), xy=(self.xmin, ypos), xytext=((self.xmax+self.xmin)*0.5, ypos),
        xycoords="data", fontsize=28, color="white") )

        self.annotations.append(
        self.ax.annotate(s="",xy=(self.xmin, ypos), xytext=(self.xmax, ypos),
        arrowprops={"arrowstyle":"<->", "edgecolor":"white"} )
        )

        self.annotations.append(
        self.ax.annotate(s="%d um"%(dy*1000.0), xy=(xpos, self.ymin), xytext=(xpos, (self.ymax+self.ymin)*0.5),
        xycoords="data", fontsize=28, color="white")
        )

        self.annotations.append(
        self.ax.annotate(s="",xy=(xpos, self.ymin), xytext=(xpos, self.ymax),
        arrowprops={"arrowstyle":"<->", "edgecolor":"white"} )
        )


    def getCenter( self, R ):
        if ( self.currentCurvature == "concave" ):
            self.yc = self.y0 - R/np.sqrt( 1 + self.tanAlpha**2 )
            self.xc = self.x0 + (self.y0-self.yc)*self.tanAlpha
        else:
            self.yc = self.y0 + R/np.sqrt( 1 + self.tanAlpha**2 )
            self.xc = self.x0 + (self.y0-self.yc)*self.tanAlpha
        return (self.xc,self.yc)

    def slope( self, x1, y1 ):
        return -(x1-self.xc)/(y1-self.yc)

    def getY( self, R, x ):
        if ( self.currentCurvature == "concave" ):
            return self.yc + np.sqrt( R*R - (x-self.xc)**2 )
        return self.yc - np.sqrt( R*R - (x-self.xc)**2 )

    def addArc( self, R, angle, curvature ):
        self.currentCurvature = curvature
        # Compute center
        self.getCenter( R )

        # Set slope to the end point
        if ( curvature == "concave" ):
            self.tanAlpha = -np.tan( -np.arctan(self.tanAlpha) + angle*np.pi/180.0 )
            y1 = self.yc + R/np.sqrt( 1.0 + self.tanAlpha**2 )
        else:
            self.tanAlpha = -np.tan( - np.arctan(self.tanAlpha) - angle*np.pi/180.0 )
            y1 = self.yc - R/np.sqrt( 1.0 + self.tanAlpha**2 )
        x1 = self.xc - (y1-self.yc)*self.tanAlpha

        x = np.linspace( self.x0, x1, 100 )
        y = self.getY( R, x )
        self.ax.plot( x, y, lw=15, color=self.colors[self.current%len(self.colors)] )

        # Update slope

        # Set new start points
        self.x0 = x1
        self.y0 = y1

        if ( np.min(y) < self.ymin ):
            self.ymin = np.min(y)
        if ( np.max(y) > self.ymax ):
            self.ymax = np.max(y)

        if ( np.min(x) < self.xmin ):
            self.xmin = np.min(x)
        if ( np.max(x) > self.xmax ):
            self.xmax = np.max(x)

        if ( self.isFirst ):
            self.xmin = np.min(x)
            self.xmax = np.max(x)
            self.ymin = np.min(y)
            self.ymax = np.max(y)
            self.isFirst = False
        self.current += 1

    def show( self ):
        for an in self.annotations:
            an.remove()
        self.annotations = []
        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_ylim( self.ymin, self.ymax )
        self.customizePlot()
        self.fig.show()
        plt.show()
        #plt.show( block=False )
