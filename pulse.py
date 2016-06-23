import numpy as np
import matplotlib
import mplLaTeX as mltx
matplotlib.rcParams.update(mltx.params)
from matplotlib import pyplot as plt
import json

def main():
    ofname = "Figures/pulse20deg.pdf"
    fname = "dataPlane/MultInc20s/WithEps/realField.json"
    fnameBkg = "dataPlane/MultInc20s/bkg/realField.json"
    infile = open( fname, 'r' )
    data = json.load(infile)
    infile.close()

    infile = open( fnameBkg, 'r' )
    dataBkg = json.load(infile)
    infile.close()

    t = np.array( data["time"] )
    reflReal = np.array( data["reflected"]["real"] )
    transReal = np.array( data["transmitted"]["real"] )

    reflBkgReal = np.array( dataBkg["reflected"]["real"] )

    fig = plt.figure()
    ax = fig.add_subplot(111)
    reflectedWave = reflReal - reflBkgReal
    minBkg = np.min( reflBkgReal )
    maxReflected = np.max( reflectedWave )
    minReflected = np.min( reflectedWave )
    maxTrans = np.max( transReal ) 

    # Compute shifts
    shiftBkg = 1.05*( np.abs(maxReflected) + np.abs( minBkg ) )
    shiftTrans = 1.05*( np.abs( minReflected ) + np.abs( maxTrans ) )

    # Compute end point
    positionOfMax = np.array( [np.argmax( reflBkgReal ), np.argmax( reflectedWave ), np.argmax( transReal )] )
    endPoint = 2*np.max( positionOfMax )
    ax.plot( t[:endPoint], reflectedWave[:endPoint], color='black')
    ax.plot( t[:endPoint], reflBkgReal[:endPoint]+shiftBkg, color='black')
    ax.plot( t[:endPoint], transReal[:endPoint]-shiftTrans, color='black')

    x1, x2 = ax.get_xlim()
    xpos = 0.75*(x2-x1)+x1
    ax.text( xpos, shiftBkg + 0.2*maxReflected, "Incident" )
    ax.text( xpos, 0.2*maxReflected, "Reflected")
    ax.text( xpos, -shiftTrans+0.2*maxReflected, "Transmitted")
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "Amplitude" )
    fig.savefig( ofname, bbox_inches="tight")

if __name__ == '__main__':
    main()
