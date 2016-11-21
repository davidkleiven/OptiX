import circularWG as cw
import sys
import numpy as np

def main(argv):
    wg = cw.CurvedWaveguide()
    wg.R = 40E6
    wg.width = 1200.0
    wg.kappa = 1.0

    wg.solve(100, 1.8*np.pi/180.0 )
    wg.plot()
    wg.rayDensity()

if __name__ == "__main__":
    main(sys.argv[1:])
