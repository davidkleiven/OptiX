from scipy import io
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

def main():
    fname = "data/delta_beta_in_sphere_updated.mat"
    matfile = io.loadmat(fname)
    beta = np.array( matfile["beta_afo_r"][0] )
    delta = np.array( matfile["delta_afo_r"][0] )
    r = np.array( matfile["rv_for_delta"][0] )*1E9
    fname = "data/radialSphereDist.csv"
    np.savetxt("data/radialSphereDist.csv", np.vstack([r,delta,beta]).T, header="r (nm), delta, beta", delimiter=",")
    print ("CSV file written to %s"%(fname))


    # Interpolate data to regular grid
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(r,delta, color="red")
    ax2 = ax.twinx()
    ax2.plot(r,beta, color="blue")
    plt.show()

if __name__ == "__main__":
    main()
