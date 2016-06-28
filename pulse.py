import numpy as np
import matplotlib
import mplLaTeX as mltx
matplotlib.rcParams.update(mltx.params)
from matplotlib import pyplot as plt
import json
import sys

'''
Example call:
python pulse.py dataPlane/MultInc20s Figures

This call looks for background files in:
dataPlane/MultInc20s/bkg

file with dielectric slab in
dataPlane/MultInc20s/bkg

The figure will be saved in 
Figures
'''
def main(argv):
    if ( len(argv) != 2 ):
        print "[pulse] Usage: python pulse.py <ddir> <figdir>"
        return
    ddir = argv[0]
    fdir = argv[1]
    positionOfFirstDigit = -1
    digitIsFound = False
    positionOfPolarisation = -1
    for i in range(0, len(ddir)):
        if ( ddir[i].isdigit() and not digitIsFound ):
            positionOfFirstDigit = i
            digitIsFound = True
        if ( (not ddir[i].isdigit()) and digitIsFound ):
            positionOfPolarisation = i
            break
    if ( positionOfPolarisation == -1 or positionOfFirstDigit == -1 ):
        ofname = fdir + "/pulseUnknownAngleAndPolarisation.pdf"
    else:
        ofname = fdir + "/pulse%s.pdf"%(ddir[positionOfFirstDigit:positionOfPolarisation+1])
            
    fname = ddir+"/WithEps/realField.json"
    fnameBkg = ddir+"/bkg/realField.json"
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
    print "[pulse] Figure written to %s"%(ofname)

if __name__ == '__main__':
    main(sys.argv[1:])
