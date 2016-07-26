USE_LATEX = False
import numpy as np
import matplotlib
if ( USE_LATEX ):
    import mplLaTeX as mltx
    matplotlib.rcParams.update(mltx.params)
from matplotlib import pyplot as plt
from scipy import signal
import json
import sys

HIGH_ANGLE = "dataPlane/GeometryClass85s/bkg/realField.json"
LOW_ANGLE = "dataPlane/GeometryClass5s/bkg/realField.json"

def envelope( sig ):
    hTransform = signal.hilbert( sig )
    return np.abs( hTransform )

def main():
    try:
        infile = open( HIGH_ANGLE, 'r' )
    except:
        print ("Could not open file %s"%(LOW_ANGLE) )
        return 
    
    highData = json.load(infile)
    infile.close()

    try:
        infile = open( LOW_ANGLE, 'r' )
    except:
        print ("Could not open file %s"%(HIGH_ANGLE))
        return
    lowData = json.load(infile)
    infile.close()

    # Plot pulses
    highAngleSignal = np.array( highData["transmitted"]["real"] )
    highEnv = envelope( highAngleSignal )
    lowAngleSignal = np.array( lowData["transmitted"]["real"] )
    lowEnv = envelope( lowAngleSignal )
    highSkewness = np.mean( np.abs(highAngleSignal)**3 )/np.std( highAngleSignal )**3
    lowSkewness = np.mean( np.abs( lowAngleSignal )**3 )/np.std( lowAngleSignal )**3

    # Low angle signal envelope
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot( lowData["time"], np.array( lowData["transmitted"]["real"] ), 'k' )
    ax.plot( lowData["time"], lowEnv, 'b' )
    
    # High angle signal and envelope
    fig2 = plt.figure() 
    ax2 = fig2.add_subplot(111)
    ax2.plot( highData["time"], highAngleSignal, 'k' )
    ax2.plot( highData["time"], highEnv, 'b')
    
    # Plot envelopes and mirrors
    lowEnv /= np.max(lowEnv)
    mirrorLow = lowEnv[::-1]
    posMaxForward = np.argmax( lowEnv )
    posMaxMirror = np.argmax( mirrorLow )  
    mirrorLow = np.roll( mirrorLow, posMaxForward-posMaxMirror-1 )
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot( mirrorLow - lowEnv)
    
    highEnv /= np.max(highEnv)
    mirrorHigh = highEnv[::-1]
    posMaxForward = np.argmax( highEnv )
    posMaxMirror = np.argmax( mirrorHigh )
    mirrorHigh = np.roll( mirrorHigh, posMaxForward-posMaxMirror-1 )
    fig4 = plt.figure() 
    ax4 = fig4.add_subplot(111)
    ax4.plot( highEnv-mirrorHigh )
    if ( not USE_LATEX ):
        plt.show()

if __name__ == '__main__':
    main()
