import numpy as np
from matplotlib import pyplot as plt
import os
import sys
from xrt.backends.raycing import materials as rm

def main():
    film = rm.Material("Si", kind="thin mirror", rho=2.32, t=100E-9)
    E = 7000
    print ("Refractive index:")
    print film.get_refractive_index(E)

if __name__ == "__main__":
    main()
