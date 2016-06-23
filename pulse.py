import numpy as np
from matplotlib import pyplot as plt
import json

def main():
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
    endPoint = 2.0*np.max( positionOfMax )
    ax.plot( t[:endPoint], reflectedWave[:endPoint], color='black')
    ax.plot( t[:endPoint], reflBkgReal[:endPoint]+shiftBkg, color='black')
    ax.plot( t[:endPoint], transReal[:endPoint]-shiftTrans, color='black')
    ax.set_xlabel( "Time (s)" )
    ax.set_ylabel( "Amplitude" )
    fig.show()
    plt.show()

if __name__ == '__main__':
    main()
