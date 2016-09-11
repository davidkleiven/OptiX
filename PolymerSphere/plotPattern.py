import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colors
import json

def formFactor( n, qR ):
    return (n**2 - 1.0)*(np.sin(qR) - qR*np.cos(qR))/(qR**3)

def scatteringPattern( n, qR ):
    return np.abs(1.0+formFactor(n,qR))**2

def main():
    filename = "data/overview0.json"
    infile = open(filename, 'r' )
    overview = json.load(infile)
    infile.close()

    data = np.fromfile(overview["ScatteredField"], dtype=np.float64)
    dataTot = np.fromfile(overview["TotalField"], dtype=np.float64)
    print ("Detector position (unit R): %.1f"%(overview["Detector"]["z"]))
    
    xmin = overview["Detector"]["min"]
    xmax = overview["Detector"]["max"]
    x = np.linspace(xmin,xmax, overview["Detector"]["pixels"])
    y = np.linspace(xmin,xmax, overview["Detector"]["pixels"])
    theta = np.arctan(x/overview["Detector"]["z"])
    qR = 2.0*overview["kR"]*np.sin(theta/2.0)

    X,Y = np.meshgrid(x,y)
    data = data.reshape((overview["Detector"]["pixels"],-1))
    dataTot = dataTot.reshape((overview["Detector"]["pixels"],-1))
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.contourf(X, Y, data, 200, cmap="gist_heat")

    '''
    figT = plt.figure()
    axT = figT.add_subplot(1,1,1)
    axT.contourf(X, Y, dataTot, 200, cmap="gist_heat")
    '''

    # Extract data through the center
    centerLine = data[int(overview["Detector"]["pixels"]/2)-1,:]
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    ax2.plot(x, centerLine/np.max(centerLine), 'k')

    pattern = np.abs(formFactor( 1-1E-5+1j*1E-6, qR ))**2
    ax2.plot( x, pattern/np.max(pattern) ) 
    #ax2.set_yscale('log')
    plt.show()

if __name__ == "__main__":
    main()
