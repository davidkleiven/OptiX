import sys
sys.path.append("/home/david/Documents/pymiecoated")
from pymiecoated import Mie
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colors 
import json
import h5py

HELP_MSG = "Usage: python plotPattern.py [--file=overview.json --help]\n"
HELP_MSG += "--help - Print this message\n"
HELP_MSG += "--file - overview file contianing information about other files\n"
HELP_MSG += "         if not present it will use the one in *data* folder\n"

def formFactor( n, qR ):
    return (n**2 - 1.0)*(np.sin(qR) - qR*np.cos(qR))/(qR**3)

def scatteringPattern( n, qR ):
    return np.abs(1.0+formFactor(n,qR))**2

def infilename( filenameInJson, folder ):
    # Verify that there is no folder present in front of the filename
    fname = filenameInJson.split("/")[-1]
    return folder+"/"+fname

def main(argv):
    filename = "data/overview0.json" 
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            filename = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print HELP_MSG
            return
    fnameSplit = filename.split("/")
    folder = fnameSplit[0]
    for i in range(1, len(fnameSplit)-1):
        print i
        folder += "/"
        folder += fnameSplit[i]
    infile = open(filename, 'r' )
    overview = json.load(infile)
    infile.close()
    eps = overview["eps"]["real"] + 1j*overview["eps"]["imag"]
    n = np.sqrt(eps)

    data = np.fromfile(infilename(overview["ScatteredField"],folder), dtype=np.float64)
    dataTot = np.fromfile(infilename(overview["TotalField"],folder), dtype=np.float64)
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
    with h5py.File(infilename(overview["XpointsCenter"],folder), 'r') as hf:
        print(hf.keys())
        posVec = hf["rVec"].value
    xCenter = posVec[:,0]
    with h5py.File(infilename(overview["FieldCenter"], folder), 'r') as hf:
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
    ax2.plot( thetaDeg, intensity/np.max(intensity), 'ko', fillstyle="none", label="BEM")

    pattern = np.abs(formFactor( n, qR ))**2
    ax2.plot( thetaDeg, pattern*np.cos(theta)**2/np.max(pattern), label="1st Born approx" ) 
    ax2.plot( thetaDeg, S1*np.cos(theta)**2/np.max(S1), label="$S_1$")
    ax2.plot( thetaDeg, S2*np.cos(theta)**2/np.max(S2), label="$S_2$")
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
    main(sys.argv[1:])
