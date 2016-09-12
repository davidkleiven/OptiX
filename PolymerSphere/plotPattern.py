import sys
sys.path.append("/home/david/Documents/pymiecoated")
from pymiecoated import Mie
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colors 
import json
import h5py

def formFactor( n, qR ):
    return (n**2 - 1.0)*(np.sin(qR) - qR*np.cos(qR))/(qR**3)

def scatteringPattern( n, qR ):
    return np.abs(1.0+formFactor(n,qR))**2

def main():
    filename = "data/overview0.json"
    infile = open(filename, 'r' )
    overview = json.load(infile)
    infile.close()
    eps = overview["eps"]["real"] + 1j*overview["eps"]["imag"]
    n = np.sqrt(eps)

    data = np.fromfile(overview["ScatteredField"], dtype=np.float64)
    dataTot = np.fromfile(overview["TotalField"], dtype=np.float64)
    print ("Detector position (unit R): %.1f"%(overview["Detector"]["z"]))
    
    xmin = overview["Detector"]["min"]
    xmax = overview["Detector"]["max"]
    x = np.linspace(xmin,xmax, overview["Detector"]["pixels"])
    y = np.linspace(xmin,xmax, overview["Detector"]["pixels"])
    theta = np.arctan(x/overview["Detector"]["z"])
    qR = 2.0*overview["kR"]*np.sin(theta/2.0)

    # Exact solution
    mie = Mie(x=overview["kR"], eps=eps, mu=1.0)
    S1 = np.zeros(len(theta))
    S2 = np.zeros(len(theta))
    for i in range(0, len(theta)):
        S1[i] = np.abs(mie.S12(np.cos(theta[i]))[0])**2
        S2[i] = np.abs(mie.S12(np.cos(theta[i]))[1])**2

    X,Y = np.meshgrid(x,y)
    data = data.reshape((overview["Detector"]["pixels"],-1))
    dataTot = dataTot.reshape((overview["Detector"]["pixels"],-1))
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.contourf(X, Y, data, 200, cmap="gist_heat")

    figT = plt.figure()
    axT = figT.add_subplot(1,1,1)
    axT.contourf(X, Y, dataTot, 200, cmap="gist_heat")

    # Extract data through the center
    with h5py.File(overview["XpointsCenter"], 'r') as hf:
        print(hf.keys())
        posVec = hf["rVec"].value
    xCenter = posVec[:,0]
    with h5py.File(overview["FieldCenter"], 'r') as hf:
        print (hf.keys()) 
        fields = hf["Fields"].value
    
    E = fields[:3,:]
    H = fields[3:,:]
    Ec = E[:,::2] + 1j*E[:,1::2]
    Hc =  H[:,::2] + 1j*H[:,1::2]
    S = np.cross(Ec, Hc.conj(), axisa=0, axisb=0)
    intensity = np.sum(np.abs(S)**2, axis=1)
    
    thetaDeg = theta*180.0/np.pi
    centerLine = data[int(overview["Detector"]["pixels"]/2)-1,:]
    indxmax = data.argmax()
    row = indxmax%overview["Detector"]["pixels"]
    centerline = data[row,:]
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
#    ax2.plot(thetaDeg, centerLine/np.max(centerLine), 'k')
    ax2.plot( thetaDeg, intensity/np.max(intensity), 'k', label="BEM")

    pattern = np.abs(formFactor( n, qR ))**2
    #ax2.plot( thetaDeg, pattern*np.cos(theta)**2/np.max(pattern) ) 
    ax2.plot( thetaDeg, S1*np.cos(theta)**2/np.max(S1), label="$S_1$")
    #ax2.plot( thetaDeg, S2*np.cos(theta)**2/np.max(S2))
    #ax2.plot( thetaDeg, np.cos(theta)**2 )
    #ax2.set_yscale('log')
    ax2.legend(frameon=False)
    ax2.set_xlabel("Scattering angle (deg)")
    ax2.set_ylabel("Normalised scattering amplitude")
    #plotAllLines(data)
    plt.show()

def plotAllLines(data):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(0, len(data[0,:])):
        ax.plot(data[i,:], label="%d"%(i))
    ax.legend()
    fig.show()
if __name__ == "__main__":
    main()
