import numpy as np

class Border:
    def __init__(self):
        self.x1 = []
        self.z1 = []
        self.x2 = []
        self.z2 = []

class WaveGuideBorders:
    def __init__(self):
        self.borders = []
        self.lw = 0.2

    def addBorder( self, x1, z1, x2, z2 ):
        border = Border()
        border.x1 = x1
        border.z1 = z1
        border.x2 = x2
        border.z2 = z2
        self.borders.append(border)

    def visualize(self, ax):
        for border in self.borders:
            ax.plot( border.z1/1E3, border.x1, color="black", lw=self.lw )
            ax.plot( border.z2/1E3, border.x2, color="black", lw=self.lw )
