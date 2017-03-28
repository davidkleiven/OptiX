import numpy as np

def main():
    N = 50
    R = 20
    data = np.zeros((N,N,N), dtype=np.uint8)
    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                rSq = (i-N/2)**2 + (j-N/2)**2 + (k-N/2)**2
                if ( rSq < R*R):
                    data[i,j,k] = 1

    data = data.ravel()
    data.tofile( "data/sphereVoxelGeo_216_%d_%d_%d.raw"%(N,N,N) )

if __name__ == "__main__":
    main()
