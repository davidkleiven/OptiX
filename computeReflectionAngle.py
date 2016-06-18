import numpy as np
import mplLaTeX
import matplotlib
matplotlib.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
from scipy import stats
import json

FOLDERS = ["dataPlane/MultInc5/WithEps", "dataPlane/MultInc20/WithEps", "dataPlane/MultInc45/WithEps", \
"dataPlane/MultInc75/WithEps", "dataPlane/MultInc85/WithEps", "dataPlane/MultInc5p/WithEps", \
"dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", \
"dataPlane/MultInc75p/WithEps", "dataPlane/MultInc85p/WithEps"]

FOLDERS=["dataPlane/MultInc5/WithEps", "dataPlane/MultInc20/WithEps"]
def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    tanReflTimesTanAngle = None
    tanAngle = None
    for folder in FOLDERS:
        fname = folder+"/realFieldFourier.json"
        try:
            infile = open(fname,'r')
            data = json.load(infile)
            infile.close()
        except:
            print ("Could not find file %s"%(fname))
            continue
        distanceFromPlane = data["geometry"]["sourcePosition"] - data["geometry"]["slabPosition"]
        phase = np.array( data["reflection"]["phase"] )
        tanAngleNew =  np.tan( np.array( data["angle"] )*np.pi/180.0) 
        k = 2.0*np.pi*np.array( data["frequency"] )
        tanReflTimesTanAngleNew = phase*np.sqrt(1.0+tanAngleNew**2)/(k*distanceFromPlane) - tanAngleNew*tanAngleNew 

        if ( tanReflTimesTanAngle is None ):
            tanReflTimesTanAngle = np.zeros(len(tanReflTimesTanAngleNew))
            tanAngle = np.zeros(len(tanAngleNew))
            tanReflTimesTanAngle = tanReflTimesTanAngleNew
            tanAngle = tanAngleNew
        else:
            tanReflTimesTanAngle = np.append(tanReflTimesTanAngle, tanReflTimesTanAngleNew)
            tanAngle = np.append(tanAngle, tanAngleNew)
    #slope, interscept, pvalue, rvalye, stderr = stats.linregress(np.log(tanAngle), np.log(tanReflTimesTanAngle))
    ax.plot( tanAngle, tanReflTimesTanAngle, '.', color="black", markersize=1, fillstyle="none")
    #ax.plot( tanAngle, np.exp(interscept)*tanAngle**slope, color='black')
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel("$\\tan \\theta_i$")
    ax.set_ylabel("$\\tan \\theta_i \\tan \\theta_r$")
    #print ("Exponent: %.2f"%(slope))
    fig.savefig("Figures/tanReflection.pdf", bbox_inches="tight")

if __name__ == "__main__":
    main()
