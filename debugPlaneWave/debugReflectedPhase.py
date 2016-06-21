import numpy as np
from matplotlib import pyplot as plt
import json

FNAME_RUN = "dataPlane/MultInc20s/WithEps/realField.json"
FNAME_BKG = "dataPlane/MultInc20s/bkg/realField.json"
ANGLE = 20.0
def main(): 
    infile = open(FNAME_RUN, 'r')
    dataRun = json.load(infile)
    infile.close()

    infile = open(FNAME_BKG, 'r')
    dataBkg = json.load(infile)
    infile.close()

    bkg = np.array( dataBkg["reflected"]["real"] )
    run = np.array( dataRun["reflected"]["real"] )
    run -= bkg
    # Compute FFT of the run signal
    ftRun = np.fft.rfft(run)
    freq = np.fft.rfftfreq( len(run), d=(dataRun["time"][1] - dataRun["time"][0]))
    
    # Compute FFT of the bkg signal
    ftBkg = np.fft.rfft(bkg)

    r = ftRun/ftBkg
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot( freq, np.angle(r) )
    ax.plot( freq, np.angle(ftBkg), '.')
    ax.plot( freq, np.angle(ftRun), 'x')
    fig.show()

    figAbs = plt.figure()
    axAbs = figAbs.add_subplot(111)
    axAbs.plot( freq, np.abs(ftBkg)**2 )
    axAbs.plot( freq, np.abs(ftRun)**2 )
    figAbs.show()

    figDiff = plt.figure()
    axDiff = figDiff.add_subplot(111)
    fc = 0.3
    expected = 6.45*5.0*np.sqrt( freq**2 - (fc*np.sin( ANGLE*np.pi/180.0))**2)
    expected = expected%(2.0*np.pi) - np.pi
    axDiff.plot( freq, expected )
    axDiff.plot( freq, -np.angle(r))
    figDiff.show()
    plt.show()

if __name__ == '__main__':
    main()
