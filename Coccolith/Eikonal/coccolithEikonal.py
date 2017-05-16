import numpy as np
import skfmm
from matplotlib import pyplot as plt

def main():
    fname = "../data/cocco8cv4Rotated_216_182_249_253.raw"
    data = np.fromfile( fname, dtype=np.uint8 )
    data = data.reshape((182,249,253), order="F")

    subtractPlaneAverage = True
    phi = np.ones(data.shape)
    #phi[1:phi.shape[0]-1, 1:phi.shape[1]-1, 1:phi.shape[2]-1] = -1
    phi[-1, :,:] = -1

    threshold = 80
    refractiveIndex = 10.0
    speed = np.ones(data.shape)
    speed[data>threshold] = 1.0/refractiveIndex


    travelTime = skfmm.travel_time(phi, speed, dx=21.6)

    fig = plt.figure()
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    centerX = int(travelTime.shape[0]/2)
    centerY = int(travelTime.shape[1]/2)
    centerZ = int(travelTime.shape[2]/2)

    travelTime -= travelTime.min()
    # Scale data to uint8 and save as raw
    travelTime *= 255/travelTime.max()
    travelTime = travelTime.astype(np.uint8)

    ax1.contour( travelTime[centerX,:,:] )
    ax2.contour( travelTime[:,centerY,:] )
    ax3.contour( travelTime[:,:,centerZ] )
    plt.show()

    travelTime = travelTime.ravel(order="F")

    ttimefname = fname.split(".raw")[0]+"TravelTime.raw"
    travelTime.tofile(ttimefname)
    print ("Travel time written to %s"%(ttimefname))

if __name__ == "__main__":
    main()
