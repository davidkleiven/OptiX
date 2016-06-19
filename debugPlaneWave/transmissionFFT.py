import numpy as np
from matplotlib import pyplot as plt
import json

FNAME_RUN = "dataPlane/MultInc20p/WithEps/realField.json"
FNAME_BKG = "dataPlane/MultInc20p/bkg/realField.json"

def main():
    infile = open(FNAME_RUN, 'r')
    dataRun = json.load(infile)
    infile.close()

    infile = open(FNAME_BKG, 'r')
    dataBkg = json.load(infile)
    infile.close()

    # Compute FFT of the run signal
    ftRun = np.fft.rfft(np.array(dataRun["transmitted"]["real"]))
    
    # Compute FFT of the bkg signal
    ftBkg = np.fft.rfft(np.array(dataBkg["transmitted"]["real"]))
    
    assert len(ftBkg) == len(ftRun)
    
    figFT = plt.figure()
    axFT = figFT.add_subplot(111)
    axFT.plot( abs(ftRun), color='red', lw=4, label="Run")
    axFT.plot( abs(ftBkg), color='blue', label="Bkg")
    axFT.legend(loc="upper right", frameon=False)
    figFT.show()

    # Relative
    transmission = ftRun/ftBkg
    figT = plt.figure()
    axT = figT.add_subplot(111)
    axT.plot( abs( transmission ) )
    figT.show()
    plt.show()

if __name__ == '__main__':
    main()
    
