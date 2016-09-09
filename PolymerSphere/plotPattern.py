import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colors

def formFactor( n, kR, ):
def main():
    filename = "data/intensity_det20_n200.bin"
    data = np.fromfile(filename, dtype=np.float64)
    print np.max(data)
    print np.min(data)
    
    # Assume data set is square
    n = int(np.sqrt(len(data)))
    data = data.reshape((n,-1))
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.contourf(data, 200, cmap="hot")

    # Extract data through the center
    centerLine = data[int(n/2),:]
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    ax2.plot(centerLine, 'k')
    
    plt.show()

if __name__ == "__main__":
    main()
