import numpy as np
import matplotlib as mpl
#mpl.rcParams.update( mp.params )
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams["font.size"] = 28
mpl.rcParams["axes.linewidth"] = 0.1
import subprocess
import eigenModes as eig
from matplotlib import pyplot as plt

class StraightEigenmodes(eig.Eigenmodes):
    def __init__(self):
        eig.Eigenmodes.__init__(self)
        self.potStrength = 8E-5

    def curvatureInducedPotential( self, x ):
        return -self.potStrength*x

    def computePerturbingWavefunction( self, mode, nModes ):
        data = np.zeros( len(self.modes[0].profile ))
        for i in range(0, nModes ):
            if ( i == mode ):
                continue
            x = np.linspace(self.modes[i].xmin, self.modes[i].xmax, len(data))
            values = self.modes[i].profile*self.modes[mode].profile*self.curvatureInducedPotential(x)
            integral = np.trapz( values )
            norm1 = np.trapz( self.modes[i].profile**2 )
            norm2 = np.trapz( self.modes[mode].profile**2 )
            integral /= np.sqrt( norm1*norm2 )
            data += integral*self.modes[i].profile/(self.modes[mode].eigenvalue-self.modes[i].eigenvalue)
        return data

    def addShade( self, ax, wgstart, wgend ):
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        width1 = wgstart-xmin
        height = ymax-ymin
        width2 = xmax-wgend
        color = "#d3d3d3"
        R1 = mpl.patches.Rectangle((xmin,ymin), width1, height, facecolor=color, edgecolor="none")
        R2 = mpl.patches.Rectangle((wgend,ymin), width2, height, facecolor=color, edgecolor="none")
        ax.add_patch(R1)
        ax.add_patch(R2)
        return ax

    def plotCorrection3Lowest( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        #ax2 = ax.twinx()
        colors = ["#e41a1c", "#377eb8", "#4daf4a"]
        for i in range(0, 2 ):
            xmin = self.modes[i].xmin
            xmax = self.modes[i].xmax
            x = np.linspace( xmin, xmax, len(self.modes[i].profile))
            correction = self.computePerturbingWavefunction( i, 13 )
            ax.plot( x, correction, color=colors[i], label="%d"%(i+1))
            ax.plot( x, correction+self.modes[i].profile, color=colors[i], lw=4, ls="--")
        self.addShade( ax, -100.0, 0.0 )
        ax.legend(loc="upper right", labelspacing=0.05, frameon=False)
        ax.set_xlabel("\$x\$ (nm)")
        ax.set_ylabel("Amplitude (a.u.)")
        fname = "Figures/firstOrderPerturb.svg"
        psname = "Figures/firstOrderPerturb.ps"
        fig.savefig(fname)
        subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname])
        print ("Figure written to %s and %s"%(fname, psname))
        plt.show()
